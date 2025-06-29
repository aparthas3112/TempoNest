#pragma once

#include <string>
#include <memory>
#include <chrono>
#include "../models/model.h"
#include "../samplers/sampler.h"

/**
 * @brief Tracks timing and generates run summary files for TempoNest analyses
 * 
 * Records start/end times, model configuration, and sampling parameters
 * Creates a human-readable summary file in the output directory
 */
class run_summary_t {
public:
    /**
     * @brief Initialize run summary tracking
     * @param output_directory Base output directory for summary file
     */
    explicit run_summary_t(const std::string& output_directory);
    
    /**
     * @brief Record the start time of the analysis
     */
    void start_timing();
    
    /**
     * @brief Record the end time and write summary file
     * @param model The model used in the analysis
     * @param sampler_settings The sampler configuration
     */
    void finish_and_write_summary(std::shared_ptr<model_t> model, 
                                  const sampler_settings_t& sampler_settings);
    
    /**
     * @brief Write parameter names file for plotting tools
     * @param model The model used in the analysis
     * @param sampler_settings The sampler configuration
     */
    void write_parameter_names(std::shared_ptr<model_t> model,
                              const sampler_settings_t& sampler_settings) const;

private:
    std::string output_dir_;
    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point end_time_;
    bool timing_started_;
    
    /**
     * @brief Generate the summary file content
     */
    std::string generate_summary_content(std::shared_ptr<model_t> model,
                                       const sampler_settings_t& sampler_settings,
                                       double runtime_seconds) const;
    
    /**
     * @brief Get current timestamp as formatted string
     */
    std::string get_timestamp() const;
    
    /**
     * @brief Format runtime duration into human-readable string
     */
    std::string format_runtime(double seconds) const;
};