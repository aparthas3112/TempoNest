#include "sampler.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include "../logger.h"
#include "multinest.h"

bool sampler_settings_t::validate() const
{
    try {
        // Get just the directory part by finding the last separator
        size_t last_sep = output_root.find_last_of("/\\");
        if (last_sep == std::string::npos) {
            // No directory structure to create
            return true;
        }

        // Extract just the directory path without the final prefix
        std::string dir_path = output_root.substr(0, last_sep);
        std::filesystem::path dir(dir_path);
        dir = dir.lexically_normal();

        // Create the directory structure if it doesn't exist
        if (!std::filesystem::exists(dir)) {
            if (!std::filesystem::create_directories(dir)) {
                logger::log_error("Failed to create directory structure: " + dir.string());
                return false;
            }
        }

        return std::filesystem::is_directory(dir);
    } catch (const std::filesystem::filesystem_error&) {
        logger::log_error("Unable to create optimiser output directory: " + output_root);
        return false;
    }
}

void sampler_settings_t::load_from_json(const json_node_t& json)
{
    if (auto dir = json.get_optional_value<string_t>("output_root")) {
        // Example: If json contains "results/newtest-" and pulsar name is "J0737-3039A"
        // this will set output_root to "results/newtest-J0737-3039A-"
        std::filesystem::path output_path(dir.value() + globals::pulsar->name + "-");
        output_root = output_path.string();
    }

    sample = json.get_optional_value<bool>("sample").value_or(true);
}

std::unique_ptr<sampler_t> sampler_factory_t::create(const json_node_t& config)
{
    /** @brief Extract sampler type from configuration */
    auto type = config.get_value<string_t>("type");

    if (type == "multinest") {
        auto settings = multinest_settings_t::from_json(config);
        if (!settings->validate()) {
            die("Invalid settings for MultiNest sampler");
        }
        return std::make_unique<multinest_sampler_t>(std::move(settings));
    }
    // Add new samplers here with additional else if statements

    die("Unknown sampler type: " + type);
}