#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include "../../json/json_node.h"
#include "../models/model.h"

struct parameter_stats_t {
    double mean;
    double stdev;
    double maximum_likelihood;
    double MAP;  // Maximum A Posteriori
};

// Forward declarations
class sampler_t;

/**
 * @brief Base class for sampler settings
 *
 * Provides common configuration options and interface for all sampling methods
 */
class sampler_settings_t {
public:

    sampler_settings_t() = default;
    virtual ~sampler_settings_t() = default;

    // Non-copyable but moveable
    sampler_settings_t(const sampler_settings_t&) = delete;
    sampler_settings_t& operator=(const sampler_settings_t&) = delete;

    /**
     * @brief Create settings from JSON configuration
     *
     * @param json JSON configuration node
     * @return std::unique_ptr<sampler_settings_t> Unique pointer to created settings
     */
    static std::unique_ptr<sampler_settings_t> from_json(const json_node_t& json);

    /**
     * @brief Validate the settings configuration
     *
     * @return bool True if settings are valid
     */
    virtual bool validate() const;

    /**
     * @brief Get the sampler type name
     *
     * @return string_t The type name of this sampler
     */
    virtual string_t get_type() const = 0;

    // Common settings

    /** @brief Directory for output chain files */
    string_t output_root = "./chains/";

    /** @brief Are we actually sampling */
    bool sample = true;

protected:

    /**
     * @brief Load settings from JSON configuration
     *
     * @param json JSON configuration node to load from
     */
    virtual void load_from_json(const json_node_t& json);
};

/**
 * @brief Factory class for creating different types of samplers
 *
 * This class provides a static factory method to create samplers based on configuration
 */
class sampler_factory_t {
public:

    /**
     * @brief Factory method implementation for creating samplers
     *
     * Creates and configures a specific sampler based on the type specified in the JSON config.
     * Validates settings before returning the sampler instance.
     *
     * @param config JSON configuration containing sampler type and settings
     * @return std::unique_ptr<sampler_t> Configured sampler instance
     * @throws runtime_error If sampler type is unknown or settings are invalid
     */
    static std::unique_ptr<sampler_t> create(const json_node_t& config);
};

class sampler_t {
private:

    /** @brief sampler_t-specific settings */
    std::unique_ptr<sampler_settings_t> settings_;

protected:

    sampler_settings_t& get_mutable_settings() { return *settings_; }

    /**
     * @brief Validate the output directory and create if necessary
     *
     * return bool True if directory is valid
     */
    bool create_output_directories();

public:

    explicit sampler_t(std::unique_ptr<sampler_settings_t> settings) : settings_(std::move(settings)) {}

    virtual ~sampler_t() = default;

    /**
     * @brief Deleted copy constructor and assignment operator, but allowing move semantics
     *
     * Samplers manage unique resources (settings_, likelihood calculations, output files)
     * and maintain internal state during sampling. Copying could lead to:
     * - Duplicate file handles and output corruption
     * - Multiple samplers processing the same likelihood calculations
     * - Inconsistent internal state between copies
     *
     * However, moving samplers is safe and useful for:
     * - Returning samplers from factory functions
     * - Storing samplers in containers
     * - Transferring ownership between components
     */
    sampler_t(const sampler_t&) = delete;
    sampler_t& operator=(const sampler_t&) = delete;

    /**
     * @brief Run the sampling algorithm
     *
     * @param model Shared pointer to the model to be sampled
     */
    virtual void run(std::shared_ptr<model_t> model) = 0;

    /**
     * @brief Output the results of sampling
     *
     */
    virtual void output_results() = 0;

    const sampler_settings_t& get_settings() const { return *settings_; }
};
