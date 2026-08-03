#ifndef QUERY_RESOLVER_H
#define QUERY_RESOLVER_H

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/**
 * @brief Représente une composante résolue du CDX index.
 *
 * @note Le std::string_view 'name' pointe directement vers une chaîne gérée
 * par ComponentResolver. Il reste valide UNIQUEMENT tant que l'instance de
 * ComponentResolver reste en vie et n'est pas modifiée.
 */
struct ResolvedComponent {
    std::size_t cid;
    std::string_view name;
};

class ComponentResolver {
public:
    ComponentResolver() = default;

    /**
     * @brief Enregistre une composante. Les CID doivent être attribués de façon contiguë (0, 1, 2...).
     * @throws std::runtime_error en cas de doublon de nom (insensible à la casse).
     */
    void register_component(std::size_t cid, const std::string& name);

    /**
     * @brief Résout un CID numérique ou un nom de composante vers un ResolvedComponent.
     * @throws std::runtime_error si la composante est introuvable (avec suggestions).
     */
    [[nodiscard]] ResolvedComponent resolve(const std::string& query) const;

    /**
     * @brief Récupère le nom d'origine d'une composante via son CID.
     */
    [[nodiscard]] const std::string& get_name(std::size_t cid) const;

    /**
     * @brief Nombre total de composantes enregistrées.
     */
    [[nodiscard]] std::size_t size() const noexcept {
        return cid_to_name_.size();
    }

private:
    std::vector<std::string> cid_to_name_;
    std::vector<std::string> lowercase_names_; // Cache parallèle pour les suggestions
    std::unordered_map<std::string, std::size_t> lowercase_name_to_cid_;

    [[nodiscard]] static std::string to_lower(std::string str);
    [[nodiscard]] static bool is_unsigned_integer(const std::string& str);
    [[nodiscard]] static std::size_t levenshtein_distance(const std::string& s1, const std::string& s2);

    [[nodiscard]] std::vector<std::string> get_suggestions(const std::string& query, std::size_t max_suggestions = 3) const;
};

#endif // QUERY_RESOLVER_H