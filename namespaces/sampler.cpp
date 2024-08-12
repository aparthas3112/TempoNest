#include "sampler.h"
#include "../json_loader.h"

namespace sampler {

// sampler flag chooses which sampler to use, 0 = MultiNest, 1 = PolyChord
sampler_t sampler = sampler_t::MULTINEST;

// IS: flag to use importance sampling, will be supported by upcoming release of multinest, at
// which point it should be set to 1
int importance_sampling = 0;

// modal: flag to allow multinest to search for multiple modes in the data. 1 = multimodal, 0 =
// single mode
int modal = 0;

// ceff: flag to set multinest to constant efficiency mode.  Adjusts sampling to maintain the
// efficiency set by the efr parameter. This is usefull for large dimensional problems (>
// 20dim), however the accuracy of the evidence suffers if importance sampling isn't used.
int constant_efficiency = 0;

// nlive: the number of live points used by multinest.  THe more you have the more it explores
// the parameter space, but the longer it takes to sample. As a guide depending on the
// dimensionality of the problem (this is the numer of parameters sampled after discounting
// those to be marginalised over): Less than 5 dimensions: 100 Between 5 and 10 dimensions: 200
// Between 10 and 50: 500
// More than 50: 1000
int live_points = 500;

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
double efficiency = 0.1;

// do sampling? 0 = No, 1 = Yes
bool sample = true;

// updInt: How often the output files are updated
int update_interval = 2000;

// nClsPar:  Number of parameters to cluster over when doing multi modal analysis
int num_cluster_parameters = 1;

void load_sampler(const string_t& filename)
{
    json_loader::load_sampler(filename);
    print();
}

void print()
{
    std::cout << "Sampler globals:" << std::endl;
    std::cout << "Sampler: " << (sampler == sampler_t::MULTINEST ? "MultiNest" : "PolyChord")
              << std::endl;
    std::cout << "Importance sampling: " << importance_sampling << std::endl;
    std::cout << "Modal: " << modal << std::endl;
    std::cout << "Constant efficiency: " << constant_efficiency << std::endl;
    std::cout << "Live points: " << live_points << std::endl;
    std::cout << "Efficiency: " << efficiency << std::endl;
    std::cout << "Sample: " << sample << std::endl;
    std::cout << "Update interval: " << update_interval << std::endl;
    std::cout << "Number of cluster parameters: " << num_cluster_parameters << std::endl;
}
}  // namespace sampler