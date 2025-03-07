#pragma once

#include <memory>
#include <optional>
#include "../types/basic_types.h"
#include "json/json_loader.h"

namespace globals {

// the pulsar we are evaluating
extern pulsar_t* pulsar;

extern int num_pulsars;

// global globals
extern bool debug;

// numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.
// Should only be set to 0 if all the priors are set in setTNPriors
extern int num_tempo2_its;

// useOriginalErrors - Use tempo2 errors before modification by TNEF/TNEQ
extern bool use_original_errors;

// test mode - to run unit tests
extern bool test_mode;

// json config document
extern json_loader_t config;

// are we using gpus
extern bool use_gpu;

void load_settings(const string_t& filename);

void print();
}  // namespace globals