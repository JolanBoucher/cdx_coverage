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

        void append(const std::string_view string) {
            rethrowWorkerException();

            if (string.size() > activeBuffer().size()) {
                submitActiveBuffer();
                waitUntilWriterIdle();

                writeDirectSynchronously(
                    string.data(),
                    string.size()
                );

                return;
            }

            if (active_offset_ + string.size() > activeBuffer().size()) {
                submitActiveBuffer();
            }

            std::memcpy(
                activeBuffer().data() + active_offset_,
                string.data(),
                string.size()
            );

            active_offset_ += string.size();
        }

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

            // 1. Préfixe du composant (nom + '\t')
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

            // 3. Tabulation
            buffer_begin[active_offset_++] = '\t';

            // 4. Couverture
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

            // 5. Saut de ligne
            buffer_begin[active_offset_++] = '\n';
        }

        /**
         * @brief Soumet le reliquat du buffer actif et attend l'écriture complète.
         */
        void flush() {
            rethrowWorkerException();

            submitActiveBuffer();
            waitUntilWriterIdle();

            rethrowWorkerException();
        }

        /**
         * @brief Vide tous les buffers et arrête proprement le worker d'arrière-plan.
         */
        void finish() {
            if (finished_) {
                rethrowWorkerException();
                return;
            }

            std::exception_ptr pending_exception;

            try {
                flush();
            } catch (...) {
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
        [[nodiscard]] static std::size_t validateBufferSize(const std::size_t buffer_size) {
            if (buffer_size < 128) {
                throw std::invalid_argument("TSV writer buffer must contain at least 128 bytes.");
            }
            return buffer_size;
        }

        [[nodiscard]] std::vector<char> &activeBuffer() noexcept {
            return buffers_[active_index_];
        }

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

        void writerLoop() noexcept {
            try {
                while (true) {
                    std::size_t buffer_index = 0;
                    std::size_t buffer_size = 0;

                    {
                        std::unique_lock lock(mutex_);
                        condition_.wait(lock, [this] {
                                return write_pending_ || stop_requested_;
                            }
                        );

                        if (!write_pending_ && stop_requested_) {
                            break;
                        }

                        buffer_index = pending_index_;
                        buffer_size = pending_size_;
                    }

                    stream_.write(buffers_[buffer_index].data(), static_cast<std::streamsize>(buffer_size));
                    if (!stream_) {
                        throw std::runtime_error("I/O error while writing TSV data.");
                    }

                    {
                        std::lock_guard lock(mutex_);
                        write_pending_ = false;
                        pending_size_ = 0;
                    }

                    condition_.notify_all();
                }
            } catch (...) {
                {
                    std::lock_guard lock(mutex_);
                    worker_exception_ = std::current_exception();
                    write_pending_ = false;
                    pending_size_ = 0;
                }
                condition_.notify_all();
            }
        }


        void writeDirectSynchronously(const char *data, std::size_t size) const {
            stream_.write(data, static_cast<std::streamsize>(size));
            if (!stream_) {
                throw std::runtime_error("I/O error while writing TSV data.");
            }
        }

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

    void writeCoverageTsvQuery(
        const std::filesystem::path &output_tsv,
        const std::vector<cdx::Coverage> &bp_cov_table,
        const std::string &component_name
    ) {
        // Pre-writing validation
        if (component_name.empty()) {
            throw std::invalid_argument("Component name cannot be empty.");
        }
        std::ofstream tsv(output_tsv, std::ios::binary);
        if (!tsv) {
            throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
        }

        DoubleBufferedTsvWriter writer(tsv);
        writer.append("component_name\tposition\tcoverage\n");

        const std::string component_prefix = component_name + '\t';

        for (std::size_t pos = 0; pos < bp_cov_table.size(); ++pos) {
            const cdx::Coverage coverage = bp_cov_table[pos];
            if (coverage >= cfg::NOT_IN_QUERY) continue;

            writer.appendRow(component_prefix, pos, coverage);
        }

        writer.finish();

        tsv.flush();
        if (!tsv) {
            throw std::runtime_error("Unable to flush TSV file: " + output_tsv.string());
        }
    }

    void writeCoverageTsvGlobal(
        const std::filesystem::path &output_tsv,
        const std::vector<cdx::Coverage> &flat_bp_cov_table,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names
    ) {
        // pre-writing validation
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
        for (std::size_t index = 1; index < bp_component_offsets.size(); ++index) {
            if (bp_component_offsets[index] < bp_component_offsets[index - 1]) {
                throw std::invalid_argument("Component offsets must be non-decreasing.");
            }
        }
        for (std::size_t component_id = 0; component_id < component_names.size(); ++component_id) {
            if (component_names[component_id].empty()) {
                throw std::invalid_argument("Component " + std::to_string(component_id) + " has an empty name.");
            }
        }

        std::ofstream tsv(output_tsv, std::ios::binary);
        if (!tsv) {
            throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
        }

        DoubleBufferedTsvWriter writer(tsv);

        writer.append("component_name\tposition\tcoverage\n");
        for (std::size_t component_id = 0; component_id < component_names.size(); ++component_id) {
            const std::string component_prefix = component_names[component_id] + '\t';
            const std::size_t component_start = bp_component_offsets[component_id];
            const std::size_t component_end = bp_component_offsets[component_id + 1];

            for (std::size_t flat_position = component_start; flat_position < component_end; ++flat_position) {
                const cdx::Coverage coverage = flat_bp_cov_table[flat_position];

                if (coverage >= cfg::NOT_IN_QUERY) {
                    continue;
                }
                writer.appendRow(component_prefix, flat_position - component_start, coverage);
            }
        }

        writer.finish();

        tsv.flush();
        if (!tsv) {throw std::runtime_error("Unable to flush TSV file: " + output_tsv.string());
        }
    }

} // namespace output