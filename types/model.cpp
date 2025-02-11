#include "model.h"

namespace model {

model_space_t model_space;  // Default constructed

void load_model(const string_t& filename)
{
    model_space = model_space_t(filename);
}

}  // namespace model