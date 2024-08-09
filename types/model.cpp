#include "model.h"
#include "../json_loader.h"

namespace model {

optional_element_t pl_red_noise;
optional_element_t pl_dm_noise;
optional_element_t efac;
optional_element_t equad;

void load_from_json(const string_t& filename)
{
    json_loader::load_from_json(filename);
}

}  // namespace model