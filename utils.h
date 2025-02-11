#pragma once
#include "types/basic_types.h"

/**
 * @brief Builds a JSON path by combining a base path with a key
 *
 * @param base The base path (can be empty)
 * @param key The key to append
 * @return string_t New path with key appended to base (with dot separator if base not empty)
 */
string_t make_json_path(const string_t& base, const string_t& key);

[[noreturn]] void die(const string_t& message);
