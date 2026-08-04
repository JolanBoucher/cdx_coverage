// Note: performance is now primarily limited by file I/O.
// Further code optimizations are unlikely to provide significant gains.

#include "output_coverage.h"
#include "config.h"

#include <array>
#include <charconv>
#include <condition_variable>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

namespace output {
    /**
     * @brief High-performance double-buffered TSV writer for overlapping CPU
     * formatting and background disk I/O.
     */
    class DoubleBufferedTsvWriter {
    public:
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 4 * 1024 * 1024;

        explicit DoubleBufferedTsvWriter(
            std::ofstream &stream,
            const std::size_t buffer_size = DEFAULT_BUFFER_SIZE
        )
            : stream_(stream),
              buffers_{
                  std::vector<char>(validateBufferSize(buffer_size)),
                  std::vector<char>(buffer_size)
              } {
            // Le thread worker est lancé une fois les buffers entièrement initialisés.
            worker_ = std::thread(
                &DoubleBufferedTsvWriter::writerLoop,
                this
            );
        }

        DoubleBufferedTsvWriter(const DoubleBufferedTsvWriter &) = delete;

        DoubleBufferedTsvWriter &operator=(const DoubleBufferedTsvWriter &) = delete;

        ~DoubleBufferedTsvWriter() noexcept {
            try {
                finish();
            } catch (...) {
                // Les destructeurs ne doivent jamais lever d'exceptions.
            }
        }

        /**
         * @brief Appends a string view to the active buffer, flushing or bypassing asynchronously if needed.
         */
        void append(const std::string_view string) {
            rethrowWorkerException();

            if (string.size() > activeBuffer().size()) {
                submitActiveBuffer();
                waitUntilWriterIdle();

                writeDirectSynchronously(string.data(), string.size());

                return;
            }

            if (active_offset_ + string.size() > activeBuffer().size()) {
                submitActiveBuffer();
            }

            std::memcpy(activeBuffer().data() + active_offset_, string.data(), string.size());

            active_offset_ += string.size();
        }

        /**
         * @brief Appends a formatted TSV data row (component, position, and coverage) to the active buffer.
         */
        void appendRow(
            const std::string_view component_prefix,
            const std::size_t position,
            const cdx::Coverage coverage
        ) {
            constexpr std::size_t MAX_POSITION_CHARS =
                    std::numeric_limits<std::size_t>::digits10 + 1;

            constexpr std::size_t MAX_COVERAGE_CHARS =
                    std::numeric_limits<cdx::Coverage>::digits10 + 1;

            const std::size_t maximum_row_size =
                    component_prefix.size() +
                    MAX_POSITION_CHARS + 1 + // '\t'
                    MAX_COVERAGE_CHARS + 1; // '\n'

            if (maximum_row_size > activeBuffer().size()) {
                throw std::length_error("TSV row exceeds writer buffer capacity.");
            }

            if (active_offset_ + maximum_row_size > activeBuffer().size()) {
                submitActiveBuffer();
            }

            char *const buffer_begin = activeBuffer().data();

            // 1. Component prefix (name + '\t')
            std::memcpy(
                buffer_begin + active_offset_,
                component_prefix.data(),
                component_prefix.size()
            );

            active_offset_ += component_prefix.size();

            // 2. Position
            const auto [position_end, position_error] = std::to_chars(
                buffer_begin + active_offset_,
                buffer_begin + activeBuffer().size(),
                position
            );

            if (position_error != std::errc{}) {
                throw std::runtime_error("Unable to format TSV position.");
            }

            active_offset_ = static_cast<std::size_t>(
                position_end - buffer_begin
            );

            // 3. Tab
            buffer_begin[active_offset_++] = '\t';

            // 4. Couverage
            const auto [coverage_end, coverage_error] = std::to_chars(
                buffer_begin + active_offset_,
                buffer_begin + activeBuffer().size(),
                coverage
            );

            if (coverage_error != std::errc{}) {
                throw std::runtime_error("Unable to format TSV coverage.");
            }

            active_offset_ = static_cast<std::size_t>(
                coverage_end - buffer_begin
            );

            // 5. add the endline char
            buffer_begin[active_offset_++] = '\n';
        }

        /** @brief Flushes the remaining active buffer contents and waits for write completion. */
        void flush() {
            rethrowWorkerException();

            submitActiveBuffer();
            waitUntilWriterIdle();

            rethrowWorkerException();
        }

        /* @brief Empty all the buffer and stop the background worker */
        void finish() {
            if (finished_) {
                rethrowWorkerException();
                return;
            }

            std::exception_ptr pending_exception;

            try { flush(); } catch (...) {
                pending_exception = std::current_exception();
            }

            {
                std::lock_guard lock(mutex_);
                stop_requested_ = true;
            }

            condition_.notify_all();

            if (worker_.joinable()) {
                worker_.join();
            }

            finished_ = true;

            if (pending_exception) {
                std::rethrow_exception(pending_exception);
            }
            rethrowWorkerException();
        }

    private:
        /** @brief Validates and ensures the TSV writer buffer size meets the minimum capacity requirement.*/
        [[nodiscard]] static std::size_t validateBufferSize(const std::size_t buffer_size) {
            if (buffer_size < 128) {
                throw std::invalid_argument("TSV writer buffer must contain at least 128 bytes.");
            }
            return buffer_size;
        }

        /* @brief Returns a reference to the currently active write buffer. */
        [[nodiscard]] std::vector<char> &activeBuffer() noexcept {
            return buffers_[active_index_];
        }

        /* @brief Submits the currently active buffer for asynchronous processing. */
        void submitActiveBuffer() {
            if (active_offset_ == 0) {
                return;
            }

            rethrowWorkerException();
            std::unique_lock lock(mutex_);
            condition_.wait(
                lock,
                [this] {
                    return !write_pending_ || worker_exception_ != nullptr;
                }
            );

            if (worker_exception_) {
                const std::exception_ptr exception = worker_exception_;
                lock.unlock();
                std::rethrow_exception(exception);
            }

            pending_index_ = active_index_;
            pending_size_ = active_offset_;
            write_pending_ = true;

            active_index_ = 1 - active_index_;
            active_offset_ = 0;

            lock.unlock();
            condition_.notify_all();
        }

        /* Blocks the calling thread until the background writer is idle (
         * no pending writes) or a worker exception occurs.
         */
        void waitUntilWriterIdle() {
            std::unique_lock lock(mutex_);
            condition_.wait(lock, [this] {
                                return !write_pending_ || worker_exception_ != nullptr;
                            }
            );

            if (worker_exception_) {
                const std::exception_ptr exception = worker_exception_;
                lock.unlock();
                std::rethrow_exception(exception);
            }
        }

        /**
 * @brief Background worker loop that asynchronously flushes pending buffers to the output stream.
 *
 * This function runs continuously in a background thread, waiting for buffers to be submitted
 * via condition variables. It writes the data to the target stream, handles I/O error checking,
 * manages thread synchronization states, and catches exceptions to safely propagate them to the main thread.
 */
        void writerLoop() noexcept {
            try {
                while (true) {
                    std::size_t buffer_index = 0;
                    std::size_t buffer_size = 0;

                    // Wait until a write is pending or a stop signal is requested
                    {
                        std::unique_lock lock(mutex_);
                        condition_.wait(lock, [this] {
                                            return write_pending_ || stop_requested_;
                                        }
                        );

                        // Exit the loop if stop was requested and no pending data remains
                        if (!write_pending_ && stop_requested_) {
                            break;
                        }

                        buffer_index = pending_index_;
                        buffer_size = pending_size_;
                    }

                    // Perform the blocking write operation to the output stream outside the lock
                    stream_.write(buffers_[buffer_index].data(), static_cast<std::streamsize>(buffer_size));
                    if (!stream_) {
                        throw std::runtime_error("I/O error while writing TSV data.");
                    }

                    // Reset pending status after a successful write
                    {
                        std::lock_guard lock(mutex_);
                        write_pending_ = false;
                        pending_size_ = 0;
                    }

                    // Notify waiting threads that the write operation has finished
                    condition_.notify_all();
                }
            } catch (...) {
                // Capture any asynchronous exceptions to rethrow them safely on the main thread
                {
                    std::lock_guard lock(mutex_);
                    worker_exception_ = std::current_exception();
                    write_pending_ = false;
                    pending_size_ = 0;
                }
                condition_.notify_all();
            }
        }


        /**
         * @brief Directly writes a data buffer to the output stream synchronously, bypassing the background queue.
         *
         * @param data Pointer to the character buffer to write.
         * @param size Number of bytes to write.
         *
         * @throws std::runtime_error If an I/O error occurs during the stream write operation.
         */
        void writeDirectSynchronously(const char *data, const std::size_t size) const {
            stream_.write(data, static_cast<std::streamsize>(size));
            if (!stream_) {
                throw std::runtime_error("I/O error while writing TSV data.");
            }
        }


        /**
         * @brief Checks for and rethrows any exception captured by the background worker thread.
         *
         * @throws std::exception_ptr Re-raises any exception that occurred during asynchronous writing.
         */
        void rethrowWorkerException() const {
            std::exception_ptr exception;
            {
                std::lock_guard lock(mutex_);
                exception = worker_exception_;
            }
            if (exception) {
                std::rethrow_exception(exception);
            }
        }


        //--- the private member of the class ---
        std::ofstream &stream_;
        std::array<std::vector<char>, 2> buffers_;

        std::size_t active_index_ = 0;
        std::size_t active_offset_ = 0;
        std::size_t pending_index_ = 0;
        std::size_t pending_size_ = 0;

        mutable std::mutex mutex_;
        std::condition_variable condition_;

        bool write_pending_ = false;
        bool stop_requested_ = false;
        bool finished_ = false;

        std::exception_ptr worker_exception_;
        std::thread worker_;
    };


    // Writes the coverage TSV query report for a specific genomic component.
    void writeCoverageTsvQuery(
        const std::filesystem::path &output_tsv,
        const std::vector<cdx::Coverage> &bp_cov_table,
        const std::string &component_name
    ) {
        // Ensure the component name is valid before opening any resources
        if (component_name.empty()) {
            throw std::invalid_argument("Component name cannot be empty.");
        }

        // Open the output file in binary mode for cross-platform newline consistency and performance
        std::ofstream tsv(output_tsv, std::ios::binary);
        if (!tsv) {
            throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
        }

        // Initialize the double-buffered writer and write the TSV header
        DoubleBufferedTsvWriter writer(tsv);
        writer.append("component_name\tposition\tcoverage\n");

        // Pre-format the component name prefix to optimize row-writing loops
        const std::string component_prefix = component_name + '\t';

        for (std::size_t pos = 0; pos < bp_cov_table.size(); ++pos) {
            const cdx::Coverage coverage = bp_cov_table[pos];

            // Skip positions that are marked as out-of-query or invalid
            if (coverage >= cfg::NOT_IN_QUERY) continue;

            writer.appendRow(component_prefix, pos, coverage);
        }

        writer.finish();

        // Ensure data is successfully written to disk and check for underlying stream errors
        tsv.flush();
        if (!tsv) {
            throw std::runtime_error("Unable to flush TSV file: " + output_tsv.string());
        }
    }

    // Write global coverage values to a TSV file using component-relative coordinates
    void writeCoverageTsvGlobal(
        const std::filesystem::path &output_tsv,
        const std::vector<cdx::Coverage> &flat_bp_cov_table,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names
    ) {
        // Validate that component boundaries describe a valid partition of
        // the flattened coverage table, with one extra boundary marking the end.
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two boundaries.");
        }
        if (component_names.size() + 1 != bp_component_offsets.size()) {
            throw std::invalid_argument("Component name count does not match component boundary count.");
        }
        if (bp_component_offsets.front() != 0) {
            throw std::invalid_argument("The first component offset must be zero.");
        }
        if (bp_component_offsets.back() != static_cast<cdx::PosBp>(flat_bp_cov_table.size())) {
            throw std::invalid_argument("The final component offset must match the flattened coverage-table size.");
        }

        // Offsets must be sorted so that each component occupies a contiguous region in flat_bp_cov_table.
        for (std::size_t index = 1; index < bp_component_offsets.size(); ++index) {
            if (bp_component_offsets[index] < bp_component_offsets[index - 1]) {
                throw std::invalid_argument("Component offsets must be non-decreasing.");
            }
        }

        // Component names are written directly to the TSV and must therefore be present for every component.
        for (std::size_t component_id = 0; component_id < component_names.size(); ++component_id) {
            if (component_names[component_id].empty()) {
                throw std::invalid_argument("Component " + std::to_string(component_id) + " has an empty name.");
            }
        }

        std::ofstream tsv(output_tsv, std::ios::binary);
        if (!tsv) {
            throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
        }

        // Buffered writer reduces the cost of emitting one TSV row per base pair.
        DoubleBufferedTsvWriter writer(tsv);

        writer.append("component_name\tposition\tcoverage\n");
        for (std::size_t component_id = 0; component_id < component_names.size(); ++component_id) {
            const std::string component_prefix = component_names[component_id] + '\t';
            const std::size_t component_start = bp_component_offsets[component_id];
            const std::size_t component_end = bp_component_offsets[component_id + 1];

            for (std::size_t flat_position = component_start; flat_position < component_end; ++flat_position) {
                const cdx::Coverage coverage = flat_bp_cov_table[flat_position];

                // Skip sentinel values that represent positions absent from the query.
                if (coverage >= cfg::NOT_IN_QUERY) {
                    continue;
                }

                // Convert from global flattened coordinates to component-relative coordinates.
                writer.appendRow(component_prefix, flat_position - component_start, coverage);
            }
        }

        writer.finish();

        tsv.flush();
        if (!tsv) {
            throw std::runtime_error("Unable to flush TSV file: " + output_tsv.string());
        }
    }
} // namespace output
