#pragma once

#include <memory>
#include "rapidjson/document.h"
#include "types/model_element.h"

/**
 * @brief A utility class for loading and parsing JSON files
 *
 * This class provides static methods for loading JSON files and parsing their contents
 * into specific data structures used in the application. It also includes utility
 * functions for safely retrieving values from JSON objects.
 */

class json_loader {
public:

    /**
     * @brief Safely get a value from a JSON object if the key is present
     *
     * This function checks if a key exists in the JSON object and retrieves its value
     * if present. It avoids using the [] operator directly on the JSON object to prevent
     * assertions when the key doesn't exist.
     *
     * @tparam T The type of the value to retrieve
     * @param json_value The JSON object to search in
     * @param key The key to look for
     * @param value The variable to store the retrieved value
     * @return true if the key exists and the value was successfully retrieved, false otherwise
     */
    template <typename T>
    static bool get_if_present(const rapidjson::Value& json_value, const char* key, T& value)
    {
        if (json_value.IsObject() && json_value.HasMember(key)) {
            return get_if_present(json_value.FindMember(key)->value, value);
        }
        return false;
    }

    /**
     * @brief Main implementation of get_if_present
     *
     * This function attempts to retrieve a value from a JSON value and store
     * it in the provided variable. It checks for type compatibility and only
     * updates the value if the types match.
     *
     * @tparam T The type of the value to retrieve
     * @param json_value The JSON value to retrieve from
     * @param value The variable to store the retrieved value
     * @return true if the value was successfully retrieved and stored, false otherwise
     */
    template <typename T>
    static bool get_if_present(const rapidjson::Value& json_value, T& value)
    {
        if constexpr (std::is_same_v<T, bool>) {
            if (json_value.IsBool()) {
                value = json_value.GetBool();
                return true;
            }
        } else if constexpr (std::is_same_v<T, int>) {
            if (json_value.IsInt()) {
                value = json_value.GetInt();
                return true;
            }
        } else if constexpr (std::is_same_v<T, unsigned>) {
            if (json_value.IsUint()) {
                value = json_value.GetUint();
                return true;
            }
        } else if constexpr (std::is_same_v<T, int64_t>) {
            if (json_value.IsInt64()) {
                value = json_value.GetInt64();
                return true;
            }
        } else if constexpr (std::is_same_v<T, uint64_t>) {
            if (json_value.IsUint64()) {
                value = json_value.GetUint64();
                return true;
            }
        } else if constexpr (std::is_same_v<T, double>) {
            if (json_value.IsNumber()) {
                value = json_value.GetDouble();
                return true;
            }
        } else if constexpr (std::is_same_v<T, float>) {
            if (json_value.IsNumber()) {
                value = static_cast<float>(json_value.GetDouble());
                return true;
            }
        } else if constexpr (std::is_same_v<T, std::string>) {
            if (json_value.IsString()) {
                value = json_value.GetString();
                return true;
            }
        } else if constexpr (std::is_same_v<T, long double>) {
            // long doubles should be saved as strings for safe loading
            if (json_value.IsString()) {
                try {
                    value = std::stold(json_value.GetString());
                    if (std::isinf(value) || std::isnan(value)) {
                        return false;
                    }
                    return true;
                } catch (const std::exception&) {
                    return false;
                }
            }
        } else {
            // This will cause a compile-time error for unsupported types
            static_assert(always_false<T>, "Unsupported type");
        }
        return false;
    }

    /**
     * @brief Parse a JSON file into a rapidjson::Document
     *
     * @param filename The name of the file to parse
     * @return rapidjson::Document The parsed JSON document
     */
    static rapidjson::Document parse_json_file(const std::string& filename);

private:

    /**
     * @brief Helper variable template for the static_assert in get_if_present
     *
     * This is used to trigger a compile-time error for unsupported types in get_if_present.
     */
    template <typename T>
    static inline constexpr bool always_false = false;
};