/**
 * @file query_resolver.h
 * @brief Component lookup and identifier-resolution utilities.
 *
 * This module provides facilities for registering graph components and
 * resolving user-supplied component identifiers. Components may be
 * referenced either by their numeric component ID (CID) or by name using
 * case-insensitive matching.
 *
 * In addition to exact resolution, the module supports suggestion
 * generation based on prefix matching and Levenshtein edit distance,
 * allowing informative error messages when a component name cannot be
 * resolved.
 *
 * The primary entry point is ComponentResolver, which maintains the
 * mapping between component IDs and component names and returns
 * ResolvedComponent objects for successful lookups.
 */

#ifndef QUERY_RESOLVER_H
#define QUERY_RESOLVER_H

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/**
 * @brief Represents a component resolved from the CDX index.
 *
 * Stores both the numeric component identifier and its associated name
 * after a successful lookup performed by ComponentResolver.
 *
 * @note The @c name string_view refers to memory owned by a
 *       ComponentResolver instance. It remains valid only while the
 *       corresponding ComponentResolver exists and its internal storage
 *       is not modified.
 */
struct ResolvedComponent {
    std::size_t cid;
    std::string_view name;
};


/**
 * @brief Registry and lookup service for graph components.
 *
 * Maintains the mapping between component IDs (CIDs) and component names,
 * providing resolution of user-supplied component identifiers. Components
 * can be looked up either by numeric CID or by name using a
 * case-insensitive comparison.
 *
 * The resolver also supports generation of similarity-based suggestions for
 * misspelled component names to improve error reporting and user experience.
 *
 * Component IDs are expected to be assigned contiguously starting at zero.
 */
class ComponentResolver {
public:
    ComponentResolver() = default;

    /**
     * @brief Register a component. Component IDs must be assigned contiguously (0, 1, 2, ...).
     * @throws std::runtime_error If a duplicate component name is detected (case-insensitive comparison).
     */
    void register_component(std::size_t cid, const std::string& name);

   /**
    * @brief Resolve a component identifier into a registered component.
    *
    * Resolution is performed in the following order:
    * 1. Numeric component ID lookup.
    * 2. Case-insensitive component name lookup.
    * 3. Suggestion generation and error reporting.
    *
    * @param query Component ID or component name supplied by the user.
    * @return The resolved component.
    *
    * @throws std::runtime_error If the component cannot be found.
    */
    [[nodiscard]] ResolvedComponent resolve(const std::string& query) const;

    /**
     * @brief Retrieve the original name associated with a component ID.
     * @param cid Component identifier.
     * @return Reference to the component's original name.
     */
    [[nodiscard]] const std::string& get_name(std::size_t cid) const;

    /**
     * @brief Return the number of registered components.
     * @return Total number of registered components.
     */
    [[nodiscard]] std::size_t size() const noexcept {
        return cid_to_name_.size();
    }

private:
    std::vector<std::string> cid_to_name_;
    std::vector<std::string> lowercase_names_; // Cache parallèle pour les suggestions
    std::unordered_map<std::string, std::size_t> lowercase_name_to_cid_;

    /* Convert a string to lowercase */
    [[nodiscard]] static std::string to_lower(std::string str);

    /* Check if an integer is unsigned*/
    [[nodiscard]] static bool is_unsigned_integer(const std::string& str);

    /**
     * @brief Compute the Levenshtein distance between two strings.
     *
     * Returns the minimum number of insertions, deletions, and substitutions
     * needed to transform one string into the other. Used for ranking
     * similarity-based component-name suggestions.
     *
     * @param s1 First string.
     * @param s2 Second string.
     * @return Edit distance between @p s1 and @p s2.
     */
    [[nodiscard]] static std::size_t levenshtein_distance(const std::string& s1, const std::string& s2);

   /**
    * @brief Generate suggested component names for an unresolved query.
    *
    * Suggestions are generated in two stages. First, component names whose
    * lowercase representation begins with the query are returned as prefix
    * matches. If insufficient matches are found, a Levenshtein-distance search
    * is performed to identify similarly spelled component names.
    *
    * Suggestions are ranked by edit distance and limited to the requested
    * maximum count.
    *
    * @param query Component name provided by the user.
    * @param max_suggestions Maximum number of suggestions to return.
    * @return Candidate component names ordered by similarity.
    */
    [[nodiscard]] std::vector<std::string> get_suggestions(const std::string& query, std::size_t max_suggestions = 3) const;
};

#endif // QUERY_RESOLVER_H