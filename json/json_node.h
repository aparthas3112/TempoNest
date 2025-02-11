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
#include "json_types.h"

// forward declare
class json_loader_t;

/**
 * @brief Wrapper class for a JSON node that provides type-safe access to its members
 *
 * This class provides a safe interface for accessing JSON object members with type checking.
 * It includes features like:
 * - Type-safe access to values with optional or throwing interfaces
 * - Support for basic types (int, double, string, bool) and arrays of these types
 * - Support for nested objects through json_node_t values
 * - Path tracking for detailed error messages
 */
class json_node_t {
public:

    explicit json_node_t(const rapidjson::Value* value, const json_loader_t& loader, string_t path) : value_(value), loader_(&loader), path_(path) {}

    /**
     * @brief Get an optional value from the node by key
     *
     * @tparam T Type to interpret the value as (int, double, string_t, bool, or json_node_t)
     * @param key Key of the value to retrieve
     * @return std::optional<T> The value if it exists and is of correct type, std::nullopt otherwise
     *
     * This method will return std::nullopt if:
     * - The key doesn't exist
     * - The node is not an object
     * If the value exists but is of the wrong type, it will call die()
     */
    template <typename T>
    std::optional<T> get_optional_value(const string_t& key) const;

    /**
     * @brief Get a required value from the node by key
     *
     * @tparam T Type to interpret the value as (int, double, string_t, bool, or json_node_t)
     * @param key Key of the value to retrieve
     * @return T The value converted to the requested type
     * @throws Will call die() if the value doesn't exist or is the wrong type
     */
    template <typename T>
    T get_value(const string_t& key) const;

    /**
     * @brief Get an optional array of values from the node by key
     *
     * @tparam T Element type for the array (int, double, string_t, bool, or json_node_t)
     * @param key Key of the array to retrieve
     * @return std::optional<std::vector<T>> Vector of values if successful, std::nullopt otherwise
     *
     * This method will return std::nullopt if:
     * - The key doesn't exist
     * - The node is not an object
     * If the value exists but is not an array or contains elements of wrong type, it will call die()
     */
    template <typename T>
    std::optional<std::vector<T>> get_optional_array(const string_t& key) const;

    /**
     * @brief Get a required array of values from the node by key
     *
     * @tparam T Element type for the array (int, double, string_t, bool, or json_node_t)
     * @param key Key of the array to retrieve
     * @return std::vector<T> Vector of values
     * @throws Will call die() if the array doesn't exist or contains invalid elements
     */
    template <typename T>
    std::vector<T> get_array(const string_t& key) const;

    /**
     * @brief Check if the node has a member with the given key
     *
     * @param key Key to check for
     * @return bool True if the node is an object and has the specified key
     */
    bool has_member(const string_t& key) const;

    /**
     * @brief Print this node's contents to the provided output stream
     *
     * @param out Output stream to print to (defaults to std::cout)
     */
    void print(std::ostream& out = std::cout) const;

    template <typename T>
    T as_value() const;

    template <typename T>
    std::vector<T> as_array() const;

private:

    friend class json_loader_t;

    template <typename T>
    static constexpr bool always_false = false;

    static void print_value(const rapidjson::Value& value, std::ostream& out);

    void mark_as_used(const string_t& path) const;

    /**
     * @brief Get node at a dotted path relative to this node
     *
     * @param inner_path The path to navigate to
     * @return std::optional<json_node_t> Node at the path if it exists, std::nullopt otherwise
     */
    std::optional<json_node_t> get_node_at_path(const string_t& inner_path) const;

    /**
     * @brief Navigate to a node using dot notation
     *
     * @param inner_path Path to navigate
     * @return const rapidjson::Value* Pointer to the value, or nullptr if not found
     */
    const rapidjson::Value* navigate_to_path(const string_t& inner_path) const;

    const rapidjson::Value* value_;
    const json_loader_t* loader_;  // parent loader for tracking access
    string_t path_;                // Path to the node from the parent loader
};

template <typename T>
T json_node_t::as_value() const
{
    return json_type_traits<T>::get(value_, path_, loader_);
}

// Template implementations for json_node_t
template <typename T>
std::optional<T> json_node_t::get_optional_value(const string_t& path) const
{
    auto node = get_node_at_path(path);
    if (!node) {
        return std::nullopt;
    }
    return node->as_value<T>();
}

template <typename T>
T json_node_t::get_value(const string_t& key) const
{
    auto value = get_optional_value<T>(key);
    if (!value.has_value()) {
        die("Required value missing or invalid: " + make_json_path(path_, key));
    }
    return value.value();
}

template <typename T>
std::vector<T> json_node_t::as_array() const
{

    if (!value_->IsArray()) {
        die("Expected array for " + path_);
    }

    std::vector<T> result;
    result.reserve(value_->Size());

    if constexpr (std::is_same_v<T, json_node_t>) {
        // For arrays of nodes, include indices in paths
        for (rapidjson::SizeType i = 0; i < value_->Size(); i++) {
            string_t element_path = path_ + "[" + std::to_string(i) + "]";
            json_node_t elem_node(&(*value_)[i], *loader_, element_path);
            result.push_back(elem_node.as_value<T>());
        }
    } else {
        // For primitive arrays, don't include indices in paths
        for (const auto& elem : value_->GetArray()) {
            json_node_t elem_node(&elem, *loader_, path_);
            result.push_back(elem_node.as_value<T>());
        }
    }

    return result;
}

template <typename T>
std::optional<std::vector<T>> json_node_t::get_optional_array(const string_t& path) const
{
    auto node = get_node_at_path(path);
    if (!node) {
        return std::nullopt;
    }
    return node->as_array<T>();
}

template <typename T>
std::vector<T> json_node_t::get_array(const string_t& key) const
{
    auto arr = get_optional_array<T>(key);
    if (!arr.has_value()) {
        die("Required array missing or invalid: " + make_json_path(path_, key));
    }
    return arr.value();
}

/**
 * @brief Type traits specialization for nested JSON objects.  This is defined in json_node.h to avoid circular dependencies.
 */
template <>
struct json_type_traits<json_node_t> {
    static bool is_valid(const rapidjson::Value* value) { return value && value->IsObject(); }

    static json_node_t get(const rapidjson::Value* value, const string_t& path, const json_loader_t* loader)
    {
        if (!is_valid(value)) {
            die("Expected object for " + path);
        }
        return json_node_t(value, *loader, path);
    }
};