#include "settings.h"
#include <logger.h>
#include "../likelihoods/gpu_functions.h"

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
    std::cout << "Settings:" << std::endl;
    std::cout << "  debug: " << debug << std::endl;
    std::cout << "  num_tempo2_its: " << num_tempo2_its << std::endl;
    std::cout << "  useOriginalErrors: " << use_original_errors << std::endl;

    if (use_gpu) {
        logger::log_info("Using GPU for calculations");
        initializeArrayFire();
    } else {
        logger::log_info("Using CPU for calculations");
    }
}
}  // namespace globals