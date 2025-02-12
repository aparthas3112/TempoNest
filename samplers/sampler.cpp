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
        std::filesystem::path dir(output_dir);
        if (!std::filesystem::exists(dir)) {
            std::filesystem::create_directories(dir);
        }
        return std::filesystem::is_directory(dir);
    } catch (const std::filesystem::filesystem_error&) {
        logger::log_error("Unable to create optimiser output directory: " + output_dir);
        return false;
    }
}

void sampler_settings_t::load_from_json(const json_node_t& json)
{
    if (auto dir = json.get_optional_value<string_t>("output_dir")) {
        output_dir = *dir;
    }

    std::string pulsarname = globals::pulsar->name;
    output_dir += "/" + pulsarname + "-";

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