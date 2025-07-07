#include "settings.h"
#include "logger.h"
#include "../likelihood/gpu_functions.h"
#include "output_formatter.h"

namespace globals {

pulsar_t* pulsar = nullptr;

int num_pulsars = 1;

bool debug = false;

int num_tempo2_its = 1;

bool use_original_errors = true;

bool test_mode = false;


#ifdef HAVE_ARRAYFIRE
bool use_gpu = true;
#else
bool use_gpu = false;
#endif

bool use_gpu_optimized = true;   // Always use optimized GPU functions

bool verbose_mode = false;  // Default to false, will be set by sampler
VerbosityLevel verbosity_level = VerbosityLevel::QUIET;  // Default to quiet

json_loader_t config;

void load_settings(const string_t& filename)
{

    config.load_json(filename);

    // Load model elements
    std::optional<json_node_t> globals = config.get_optional_value<json_node_t>("globals");
    if (globals.has_value()) {
        const auto& json_globals = globals.value();
        // Parse sampler globals as needed

        debug = json_globals.get_optional_value<bool>("debug").value_or(false);
        num_tempo2_its = json_globals.get_optional_value<int>("num_tempo2_its").value_or(1);
        use_original_errors = json_globals.get_optional_value<bool>("use_original_errors").value_or(true);
        test_mode = json_globals.get_optional_value<bool>("test_mode").value_or(false);

        // see if we want to override use_gpu
        bool use_gpu_override = json_globals.get_optional_value<bool>("use_gpu").value_or(use_gpu);
        if (use_gpu_override && !use_gpu) {
            die("GPU support not available");
        }
        use_gpu = use_gpu_override;
        
    }

    print();
}

void print()
{
    output_formatter::print_header("TEMPONEST CONFIGURATION");
    
    output_formatter::print_section("Global Settings");
    output_formatter::print_setting("Debug Mode", debug ? "ON" : "OFF", debug);
    output_formatter::print_setting("Tempo2 Iterations", std::to_string(num_tempo2_its));
    output_formatter::print_setting("Original Errors", use_original_errors ? "ON" : "OFF", use_original_errors);
    output_formatter::print_setting("Test Mode", test_mode ? "ON" : "OFF", test_mode);

    if (use_gpu) {
        logger::log_info("Using GPU for calculations");
        initializeArrayFire();
        // Get GPU device info after initialization
        output_formatter::print_setting("GPU Acceleration", "ON", true);
        output_formatter::print_setting("GPU Optimized", use_gpu_optimized ? "ON" : "OFF", use_gpu_optimized);
    } else {
        logger::log_info("Using CPU for calculations");
        output_formatter::print_setting("GPU Acceleration", "OFF", false);
    }
    
    
    output_formatter::print_separator();
}

void set_verbose_mode(bool verbose)
{
    verbose_mode = verbose;
}

void set_verbosity_level(VerbosityLevel level)
{
    verbosity_level = level;
    // Also update verbose_mode for backward compatibility
    verbose_mode = (level == VerbosityLevel::FULL);
}
}  // namespace globals