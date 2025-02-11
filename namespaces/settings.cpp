#include "settings.h"

namespace globals {

pulsar_t* pulsar = nullptr;

bool debug = false;

string_t root = "results/Example1-";

int num_tempo2_its = 1;

bool use_original_errors = true;

bool test_mode = false;

json_loader_t config;

void load_settings(const string_t& filename)
{

    config.load_json(filename);

    std::cout << "load settings" << std::endl;

    // Load model elements
    std::optional<json_node_t> globals = config.get_optional_value<json_node_t>("globals");
    if (globals.has_value()) {
        const auto& json_globals = globals.value();
        // Parse sampler globals as needed

        debug = json_globals.get_optional_value<bool>("debug").value_or(false);
        num_tempo2_its = json_globals.get_optional_value<int>("num_tempo2_its").value_or(1);
        root = json_globals.get_optional_value<string_t>("root").value_or("results/Example1-");
        use_original_errors = json_globals.get_optional_value<bool>("use_original_errors").value_or(true);
        test_mode = json_globals.get_optional_value<bool>("test_mode").value_or(false);
    }

    print();
}

void print()
{
    std::cout << "Settings:" << std::endl;
    std::cout << "  debug: " << debug << std::endl;
    std::cout << "  root: " << root << std::endl;
    std::cout << "  num_tempo2_its: " << num_tempo2_its << std::endl;
    std::cout << "  useOriginalErrors: " << use_original_errors << std::endl;
}
}  // namespace globals