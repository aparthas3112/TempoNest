#include "settings.h"
#include "../json_loader.h"

namespace globals {

pulsar_t* pulsar = nullptr;

bool debug = false;

string_t root = "results/Example1-";

int num_tempo2_its = 1;

bool use_original_errors = true;

void load_settings(const string_t& filename)
{
    json_loader::load_settings(filename);
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