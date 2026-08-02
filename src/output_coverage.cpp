#include "output_coverage.h"
#include "config.h"

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <charconv>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <vector>


namespace output {
    /**
 * @brief High-performance buffered TSV writer.
 *
 * Avoids per-line memory allocations and locale overhead by using a raw
 * internal character buffer and std::to_chars for numerical formatting.
 */
    class BufferedTsvWriter {
    public:
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 4 * 1024 * 1024; // 4 MiB

        explicit BufferedTsvWriter(std::ofstream &stream, const std::size_t buffer_size = DEFAULT_BUFFER_SIZE)
            : stream_(stream), buffer_(buffer_size) {
        }

        ~BufferedTsvWriter() {
            try {
                flush();
            } catch (...) {
            } // Destructors must not throw exceptions
        }

        /** @brief Appends a single character (e.g., '\t' or '\n'). */
        void append(const char c) {
            if (offset_ >= buffer_.size()) {
                flush();
            }
            buffer_[offset_++] = c;
        }

        /** @brief Appends a string slice / string_view. */
        void append(const std::string_view str) {
            if (offset_ + str.size() >= buffer_.size()) {
                flush();
                // Directly write if string size exceeds buffer capacity
                if (str.size() >= buffer_.size()) {
                    stream_.write(str.data(), static_cast<std::streamsize>(str.size()));
                    return;
                }
            }
            std::memcpy(buffer_.data() + offset_, str.data(), str.size());
            offset_ += str.size();
        }

        /** @brief Appends an integral value formatted using std::to_chars. */
        template<typename T>
        std::enable_if_t<std::is_integral_v<T> > append(T value) {
            // Maximum length of a 64-bit integer string representation is ~20 chars + sign
            if (offset_ + 32 >= buffer_.size()) {
                flush();
            }

            auto [ptr, ec] = std::to_chars(
                buffer_.data() + offset_,
                buffer_.data() + buffer_.size(),
                value
            );

            if (ec != std::errc{}) {
                throw std::runtime_error("Numerical conversion error in std::to_chars");
            }
            offset_ = static_cast<std::size_t>(ptr - buffer_.data());
        }

        /**  @brief Flushes buffered data to the underlying stream. */
        void flush() {
            if (offset_ > 0) {
                stream_.write(buffer_.data(), static_cast<std::streamsize>(offset_));
                offset_ = 0;
            }
        }

    private:
        std::ofstream &stream_;
        std::vector<char> buffer_;
        std::size_t offset_ = 0;
    };

    void writeCoverageTsvQuery(
        const std::filesystem::path &output_tsv,
        const std::vector<cdx::Coverage> &bp_cov_table,
        const std::string &component_name
    ) {
        if (component_name.empty()) {
            throw std::invalid_argument("Component name cannot be empty.");
        }

        std::ofstream tsv(output_tsv, std::ios::binary);
        if (!tsv) {
            throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
        }

        BufferedTsvWriter writer(tsv);

        // TSV Header
        writer.append("component_name\tposition\tcoverage\n");

        for (std::size_t pos = 0; pos < bp_cov_table.size(); ++pos) {
            const cdx::Coverage coverage = bp_cov_table[pos];
            if (coverage >= cfg::NOT_IN_QUERY) continue;

            writer.append(component_name);
            writer.append('\t');
            writer.append(pos);
            writer.append('\t');
            writer.append(coverage);
            writer.append('\n');
        }

        // Flush remaining data from buffer
        writer.flush();
    }

    void writeCoverageTsvGlobal(
    const std::filesystem::path& output_tsv,
    const std::vector<cdx::Coverage>& flat_bp_cov_table,
    const std::vector<cdx::PosBp>& bp_component_offsets,
    const std::vector<std::string>& component_names
) {
    // Validate boundary offsets array
    if (bp_component_offsets.size() < 2) {
        throw std::invalid_argument(
            "bp_component_offsets must contain at least two boundaries."
        );
    }

    // Validate that the number of names matches the number of components (offsets - 1)
    if (component_names.size() + 1 != bp_component_offsets.size()) {
        throw std::invalid_argument("Component name count does not match component boundary count."
        );
    }

    std::ofstream tsv(output_tsv, std::ios::binary);
    if (!tsv) {
        throw std::runtime_error("Unable to open TSV file: " + output_tsv.string());
    }

    BufferedTsvWriter writer(tsv);

    // Updated TSV Header
    writer.append("component_name\tposition\tcoverage\n");

    // Iterate over graph components
    for (std::size_t cid = 0; cid + 1 < bp_component_offsets.size(); ++cid) {
        const std::string& component_name = component_names[cid];

        const auto component_start = static_cast<std::size_t>(
            bp_component_offsets[cid]
        );

        const auto component_end = static_cast<std::size_t>(
            bp_component_offsets[cid + 1]
        );

#ifndef NDEBUG
        if (component_end < component_start || component_end > flat_bp_cov_table.size()) {
            throw std::out_of_range("Invalid component boundary.");
        }
#endif

        // Contiguous linear iteration per component
        for (std::size_t flat_bp = component_start; flat_bp < component_end; ++flat_bp) {
            const cdx::Coverage coverage = flat_bp_cov_table[flat_bp];
            if (coverage >= cfg::NOT_IN_QUERY) {
                continue;
            }

            const std::size_t local_pos = flat_bp - component_start;

            writer.append(component_name);
            writer.append('\t');

            writer.append(local_pos);
            writer.append('\t');

            writer.append(coverage);
            writer.append('\n');
        }
    }

    // Flush remaining data from buffer
    writer.flush();
}
} // namespace output
