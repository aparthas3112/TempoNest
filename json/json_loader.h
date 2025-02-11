// json_loader.hpp
#pragma once

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/prettywriter.h>
#include <iostream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include "../types/basic_types.h"
#include "../utils.h"
#include "json_node.h"

using json_node_array_t = std::vector<json_node_t, std::allocator<json_node_t>>;

/**
 * @brief JSON configuration loader that tracks key usage and provides type-safe access
 *
 * This class provides a safe interface for loading and accessing JSON configuration files.
 * It includes features like:
 * - Dot notation for accessing nested values (e.g., "server.config.timeout")
 * - Type-safe access to values with optional or throwing interfaces
 * - Tracking of accessed keys to identify unused configuration entries
 * - Support for basic types (int, double, string, bool) and arrays of these types
 * - Support for nested objects through json_node_t
 */
class json_loader_t {
public:

    /**
     * @brief Loads and parses a JSON file
     *
     * @param filename Path to the JSON file to load
     * @return bool False if file cannot be opened or contains invalid JSON
     */
    bool load_json(const string_t& filename, std::optional<std::reference_wrapper<std::unordered_set<string_t>>> included_files = std::nullopt);

    /**
     * @brief Get an optional value using dot notation path
     *
     * @tparam T Type to interpret the value as (int, double, string_t, bool, or json_node_t)
     * @param path Path to the value, can include dots for nested values (e.g., "server.config.timeout")
     * @return std::optional<T> The value if it exists and is of correct type, std::nullopt otherwise
     *
     * The key will be marked as "used" for validation purposes.
     */
    template <typename T>
    std::optional<T> get_optional_value(const string_t& path) const
    {
        std::optional<json_node_t> node = get_node_at_path(path);
        if (!node.has_value()) {
            return std::nullopt;
        }
        return node->as_value<T>();
    }

    /**
     * @brief Get a required value using dot notation path
     *
     * @tparam T Type to interpret the value as (int, double, string_t, bool, or json_node_t)
     * @param path Path to the value, can include dots for nested values (e.g., "server.config.timeout")
     * @return T The value converted to the requested type
     * @throws Will call die() if the value doesn't exist or is wrong type
     *
     * The key will be marked as "used" for validation purposes.
     */
    template <typename T>
    T get_value(const string_t& path) const
    {
        std::optional<json_node_t> node = get_node_at_path(path);
        if (!node.has_value()) {
            die("Required value missing: " + path);
        }
        // The node is already pointing at the final value, just get it directly
        return node->as_value<T>();
    }

    /**
     * @brief Get an optional array using dot notation path
     *
     * @tparam T Element type for the array (int, double, string_t, bool, or json_node_t)
     * @param path Path to the array, can include dots for nested values
     * @return std::optional<std::vector<T>> Vector of values if successful, std::nullopt otherwise
     *
     * The key will be marked as "used" for validation purposes.
     */
    template <typename T>
    std::optional<std::vector<T>> get_optional_array(const string_t& path) const
    {
        std::optional<json_node_t> node = get_node_at_path(path);
        if (!node.has_value()) {
            return std::nullopt;
        }
        return node->as_array<T>();
    }

    /**
     * @brief Get a required array using dot notation path
     *
     * @tparam T Element type for the array (int, double, string_t, bool, or json_node_t)
     * @param path Path to the array, can include dots for nested values
     * @return std::vector<T> Vector of values
     * @throws Will call die() if the array doesn't exist or contains invalid elements
     *
     * The key will be marked as "used" for validation purposes.
     */
    template <typename T>
    std::vector<T> get_array(const string_t& path) const
    {
        std::optional<json_node_t> node = get_node_at_path(path);
        if (!node) {
            die("Required array missing: " + path);
        }
        return node->as_array<T>();
    }

    /**
     * @brief Prints the entire document to the provided output stream
     *
     * @param out Output stream to print to (defaults to std::cout)
     */
    void print(std::ostream& out = std::cout) const;

    /**
     * @brief Prints a specific node to the provided output stream
     *
     * @param path Path to the node to print, can include dots for nested values
     * @param out Output stream to print to (defaults to std::cout)
     * @return bool False if the path doesn't exist
     *
     * The key will be marked as "used" for validation purposes.
     */
    bool print_node(const string_t& path, std::ostream& out = std::cout) const;

    /**
     * @brief Validates that all keys in the JSON document have been accessed
     *
     * @return bool False if there are unused keys (which are logged as errors)
     *
     * This method helps identify configuration entries that were never accessed,
     * which might indicate typos in keys or outdated configuration entries.
     */
    bool validate() const;

    /**
     * @brief Marks a key as having been accessed
     *
     * @param path The path that was accessed
     */
    void mark_as_used(const string_t& path) const;

    /**
     * @brief Return the full document
     *
     */
    json_node_t get_root_node() const { return json_node_t(&doc, *this, ""); }

private:

    rapidjson::Document doc;
    mutable std::unordered_set<string_t> used_paths;

    /**
     * @brief Get parent node for a dotted path
     *
     * @param path The path to navigate to
     * @return std::optional<json_node_t> Node at the path if it exists, std::nullopt otherwise
     */
    std::optional<json_node_t> get_node_at_path(const string_t& path) const;

    /**
     * @brief Recursively collects paths that haven't been accessed
     *
     * @param path_to_value Current path in the JSON hierarchy
     * @param value Current JSON value being examined
     * @param unused_paths Output vector for collecting unused paths
     */
    void collect_unused_paths(const string_t& path_to_value, const rapidjson::Value& value, std::vector<string_t>& unused_paths) const;

    /**
     * @brief Merges a JSON document into the current document
     *
     * @param other The document to merge in
     * @throws Will call die() on type conflicts or duplicate non-array values
     */
    void merge_document(const rapidjson::Document& other);

    /**
     * @brief Recursively merges two JSON values
     *
     * @param target Value to merge into
     * @param source Value to merge from
     * @param path Current path for error reporting
     * @throws Will call die() on type conflicts or duplicate non-array values
     */
    void merge_objects(rapidjson::Value& target, const rapidjson::Value& source, const string_t& path);

    /**
     * @brief Merges arrays, concatenating them if types match
     *
     * @param target Array to merge into
     * @param source Array to merge from
     * @param path Current path for error reporting
     * @throws Will call die() if array element types don't match
     */
    void merge_arrays(rapidjson::Value& target, const rapidjson::Value& source, const string_t& path);

    /**
     * @brief Processes any INCLUDE directives in the current document.  Will die on errors.
     *
     * @param base_path Base path for included file resolution
     * @param included_files Set of already included files to prevent cycles
     */
    void process_includes(const string_t& base_path, std::unordered_set<string_t>& included_files);
};

// helper functions for making and using json files
std::string create_temp_json_file(const std::string& content);
bool load_temp_json_file(const std::string& content, json_loader_t& loader);

template <typename T>
void set_if_exists(const json_node_t& json, const std::string& key, T& value)
{
    if (auto val = json.get_optional_value<T>(key)) {
        value = *val;
    }
}