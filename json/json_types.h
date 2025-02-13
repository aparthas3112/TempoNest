#pragma once

#include <rapidjson/document.h>
#include "../types/basic_types.h"
#include "../utils.h"

// forward declare
class json_loader_t;

/**
 * @brief Type traits system for JSON value type checking and conversion
 *
 * This template provides compile-time type information and runtime conversion
 * for JSON values. It centralizes type checking logic and handles type-specific
 * conversion rules.
 */
template <typename T>
struct json_type_traits {
    // Helper variable template for static_assert with dependent types
    template <typename U>
    static constexpr bool always_false = false;

    // Default implementation that fails static_assert for unsupported types
    static_assert(always_false<T>, "Unsupported JSON type");
};

/**
 * @brief Type traits specialization for integers
 */
template <>
struct json_type_traits<int> {
    static bool is_valid(const rapidjson::Value* value) { return value && value->IsInt(); }

    static int get(const rapidjson::Value* value, const string_t& path, const json_loader_t* = nullptr)
    {
        if (!is_valid(value)) {
            die("Expected Int for " + path);
        }
        return value->GetInt();
    }
};

/**
 * @brief Type traits specialization for doubles
 *
 * Handles both decimal numbers and integers, converting integers to doubles
 * automatically. This allows numbers like '1' and '1.0' to be read as doubles.
 */
template <>
struct json_type_traits<double> {
    static bool is_valid(const rapidjson::Value* value) { return value && (value->IsDouble() || value->IsInt()); }

    static double get(const rapidjson::Value* value, const string_t& path, const json_loader_t* = nullptr)
    {
        if (!is_valid(value)) {
            die("Expected Double for " + path);
        }
        return value->IsInt() ? static_cast<double>(value->GetInt()) : value->GetDouble();
    }
};

/**
 * @brief Type traits specialization for strings
 */
template <>
struct json_type_traits<string_t> {
    static bool is_valid(const rapidjson::Value* value) { return value && value->IsString(); }

    static string_t get(const rapidjson::Value* value, const string_t& path, const json_loader_t* = nullptr)
    {
        if (!is_valid(value)) {
            die("Expected String for " + path);
        }
        return value->GetString();
    }
};

/**
 * @brief Type traits specialization for booleans
 */
template <>
struct json_type_traits<bool> {
    static bool is_valid(const rapidjson::Value* value) { return value && value->IsBool(); }

    static bool get(const rapidjson::Value* value, const string_t& path, const json_loader_t* = nullptr)
    {
        if (!is_valid(value)) {
            die("Expected Bool for " + path);
        }
        return value->GetBool();
    }
};
