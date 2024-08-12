#pragma once

#include <memory>
#include <optional>
#include "../types/basic_types.h"

enum class sampler_t { MULTINEST, POLYCHORD };

namespace sampler {

// sampler flag chooses which sampler to use, 0 = MultiNest, 1 = PolyChord
extern sampler_t sampler;

// IS: flag to use importance sampling, will be supported by upcoming release of multinest, at
// which point it should be set to 1
extern int importance_sampling;

// modal: flag to allow multinest to search for multiple modes in the data. 1 = multimodal, 0 =
// single mode
extern int modal;

// ceff: flag to set multinest to constant efficiency mode.  Adjusts sampling to maintain the
// efficiency set by the efr parameter. This is usefull for large dimensional problems (>
// 20dim), however the accuracy of the evidence suffers if importance sampling isn't used.
extern int constant_efficiency;

// nlive: the number of live points used by multinest.  THe more you have the more it explores
// the parameter space, but the longer it takes to sample. As a guide depending on the
// dimensionality of the problem (this is the numer of parameters sampled after discounting
// those to be marginalised over): Less than 5 dimensions: 100 Between 5 and 10 dimensions: 200
// Between 10 and 50: 500
// More than 50: 1000
extern int live_points;

// efr: The sampling efficiency.  The Lower this number the more carefully multinest explores
// the parameter space, so the longer it takes. The value depends on the dimensionality and the
// goal. For parameter estimation it can be set higher than if an accurate Evidence value is
// required. As a rough guide: Less than 10 dimensions: 0.8 (parameter estimation), 0.3
// (Evidence evaluation) Between 10 and 20 dimensions: 0.3 (parameter estimation), 0.1 (Evidence
// evaluation) Between 20 and 50: 0.1 (parameter estimation), 0.05 (Evidence evaluation) More
// than 50: 0.05 (parameter estimation), 0.01 (Evidence evaluation) These numbers may well be
// adjusted with experience, for parameter estimation of a single modal "blob", these may be
// much lower than required. Also: In constant efficiency mode this must be set lower, to 0.05
// for D<50 and 0.01 D>50.
extern double efficiency;

// do sampling? 0 = No, 1 = Yes
extern bool sample;

// updInt: How often the output files are updated
extern int update_interval;

// nClsPar:  Number of parameters to cluster over when doing multi modal analysis
extern int num_cluster_parameters;

void load_sampler(const string_t& filename);

void print();
}  // namespace sampler