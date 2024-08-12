#pragma once

#include <memory>
#include <optional>
#include "../types/basic_types.h"

namespace settings {

// global settings
extern bool debug;

// Root of the results files,relative to the directory in which TempoNest is run. This will be
// followed by the pulsar name, and then the individual output file extensions.
extern string_t root;

// numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.
// Should only be set to 0 if all the priors are set in setTNPriors
extern int num_tempo2_its;

// useOriginalErrors - Use tempo2 errors before modification by TNEF/TNEQ
extern bool use_original_errors;

void load_settings(const string_t& filename);

void print();
}  // namespace settings