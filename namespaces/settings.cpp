#include "settings.h"
#include "../json_loader.h"

namespace globals {

pulsar_t* pulsar = nullptr;

bool debug = false;

string_t root = "results/Example1-";

int num_tempo2_its = 1;

bool use_original_errors = true;

bool test_mode = false;

void load_settings(const string_t& filename)
{

    rapidjson::Document doc = json_loader::parse_json_file(filename);

    std::cout << "load settings" << std::endl;

    // Load model elements
    if (doc.HasMember("globals")) {
        const auto& json_settings = doc["globals"];

        // Parse sampler globals as needed
        json_loader::get_if_present(json_settings, "debug", globals::debug);
        json_loader::get_if_present(json_settings, "num_tempo2_its", globals::num_tempo2_its);
        json_loader::get_if_present(json_settings, "root", globals::root);
        json_loader::get_if_present(json_settings, "use_original_errors",
                                    globals::use_original_errors);
        json_loader::get_if_present(json_settings, "test_mode", globals::test_mode);
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