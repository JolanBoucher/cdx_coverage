#include "query_resolver.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Convert a string to lowercase
std::string ComponentResolver::to_lower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return str;
}

// Security against signed value produced by the CLI args
bool ComponentResolver::is_unsigned_integer(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(), [](unsigned char c) {
        return std::isdigit(c);
    });
}

void ComponentResolver::register_component(std::size_t cid, const std::string& name) {
    std::string normalized = to_lower(name);

    if (lowercase_name_to_cid_.find(normalized) != lowercase_name_to_cid_.end()) {
        throw std::runtime_error("Duplicate component name detected in CDX: '" + name + "'");
    }

    if (cid >= cid_to_name_.size()) {
        cid_to_name_.resize(cid + 1);
        lowercase_names_.resize(cid + 1);
    }

    cid_to_name_[cid] = name;
    lowercase_names_[cid] = normalized;
    lowercase_name_to_cid_[normalized] = cid;
}

const std::string& ComponentResolver::get_name(std::size_t cid) const {
    if (cid >= cid_to_name_.size()) {
        throw std::out_of_range("Component CID " + std::to_string(cid) + " is out of bounds.");
    }
    return cid_to_name_[cid];
}

// Calculate the minimum number of step to go from one string to the other
std::size_t ComponentResolver::levenshtein_distance(const std::string& s1, const std::string& s2) {
    const std::size_t m = s1.size();
    const std::size_t n = s2.size();
    std::vector dp(m + 1, std::vector<std::size_t>(n + 1));

    for (std::size_t i = 0; i <= m; ++i) dp[i][0] = i;
    for (std::size_t j = 0; j <= n; ++j) dp[0][j] = j;

    for (std::size_t i = 1; i <= m; ++i) {
        for (std::size_t j = 1; j <= n; ++j) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + std::min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
            }
        }
    }
    return dp[m][n];
}

// When users give the wrong component name input generate suggesting based on the similarity with the input string
std::vector<std::string> ComponentResolver::get_suggestions(
    const std::string& query,
    std::size_t max_suggestions) const {

    std::string lower_query = to_lower(query);
    std::vector<std::string> suggestions;

    // Step 1: Prefer prefix matches because they are typically the most
    // relevant and avoid the cost of edit-distance calculations.
    for (std::size_t cid = 0; cid < lowercase_names_.size(); ++cid) {
        if (lowercase_names_[cid].rfind(lower_query, 0) == 0) { // Commences-with
            suggestions.push_back(cid_to_name_[cid]);
            if (suggestions.size() >= max_suggestions) {
                return suggestions;
            }
        }
    }

    // Step 2: Search for similarly spelled names using a dynamic edit-distance threshold based on query length.
    const std::size_t max_distance = std::max<std::size_t>(2, lower_query.size() / 4);
    std::vector<std::pair<std::size_t, std::string>> candidates;

    for (std::size_t cid = 0; cid < lowercase_names_.size(); ++cid) {
        std::size_t dist = levenshtein_distance(lower_query, lowercase_names_[cid]);
        if (dist <= max_distance) {
            candidates.emplace_back(dist, cid_to_name_[cid]);
        }
    }

    // Rank candidates by increasing edit distance.
    std::sort(candidates.begin(), candidates.end(), [](
            const auto& a, const auto& b) {
                return a.first < b.first;
    });

    for (std::size_t i = 0; i < std::min(max_suggestions, candidates.size()); ++i) {
        suggestions.push_back(candidates[i].second);
    }

    return suggestions;
}

// Resolve a component identifier into a registered component
ResolvedComponent ComponentResolver::resolve(const std::string& query) const {
    // Rule 1: Purely numeric input is interpreted as a component ID.
    if (is_unsigned_integer(query)) {
        std::size_t cid = std::stoull(query);
        if (cid < cid_to_name_.size() && !cid_to_name_[cid].empty()) {
            return ResolvedComponent{cid, cid_to_name_[cid]};
        }

        // Numeric identifier is outside the registered component range.
        std::ostringstream oss;
        oss << "Error: Component ID '" << cid << "' is out of range (total components: " << cid_to_name_.size() << ").\n";
        oss << "Use --inspect to list available components.";
        throw std::runtime_error(oss.str());
    }

    // Rule 2: Attempt a case-insensitive lookup by component name.
    std::string lower_query = to_lower(query);
    auto it = lowercase_name_to_cid_.find(lower_query);

    if (it != lowercase_name_to_cid_.end()) {
        std::size_t cid = it->second;
        return ResolvedComponent{cid, cid_to_name_[cid]};
    }

    // Rule 3: Build an informative error message and include suggested component names when close matches exist.
    std::ostringstream oss;
    oss << "Component '" << query << "' was not found.\n";

    auto suggestions = get_suggestions(query);
    if (!suggestions.empty()) {
        oss << "Did you mean:\n";
        for (const auto& sug : suggestions) {
            oss << "  - " << sug << "\n";
        }
    } else {
        oss << "\nUse --inspect to list available components.";
    }

    throw std::runtime_error(oss.str());
}