#include "utils.h"
#include "logger.h"

[[noreturn]] void die(const string_t& message)
{
    // Print the error message
    logger::log_error(message);

    // Throw a runtime_error instead of a custom exception
    throw std::runtime_error(message);
}

string_t make_json_path(const string_t& base, const string_t& key)
{
    return base.empty() ? key : base + "." + key;
}