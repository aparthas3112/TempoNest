#pragma once
#include "../../../eigen_config.h"

#ifdef HAVE_ARRAYFIRE
    #include <arrayfire.h>

// Function to convert Eigen matrix to ArrayFire array (double precision)
af::array eigenToAf(const Eigen::MatrixXd& eigenMat);

// Function to convert Eigen vector to ArrayFire array (double precision)
af::array eigenToAf(const Eigen::VectorXd& eigenVec);

namespace gpu_data {
// Static ArrayFire arrays
extern af::array total_matrix_;
extern bool initialized;

// Performance enhancement: pre-allocated working arrays
extern af::array NT_cache_;
extern af::array TNT_cache_;
extern af::array workspace_cache_;
extern int cached_size_;

// Initialize the static matrices
void initialize(const Eigen::MatrixXd& eigenTotalMatrix);

// Initialize with size pre-allocation for performance
void initialize_with_cache(const Eigen::MatrixXd& eigenTotalMatrix, int max_expected_size);

// Check if initialized
bool isInitialized();

// Clean up resources if needed
void cleanup();

// Pre-allocate working memory for repeated operations
void preallocate_workspace(int matrix_size);
}  // namespace gpu_data

#else
// Stub namespace when ArrayFire is not available
namespace gpu_data {
extern bool initialized;

void initialize(const Eigen::MatrixXd&);
bool isInitialized();
void cleanup();
int get_batch_size();
std::string get_gpu_info();
}  // namespace gpu_data
#endif

void initializeArrayFire();

// Original function (maintained for compatibility)
double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                      double freq_det, double timelike, double uniform_prior);

// Optimized version with pre-allocated memory and vectorized operations
double performAlgebraWithArrayFireGPU_optimized(const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                               double freq_det, double timelike, double uniform_prior);

// Batch processing for multiple parameter sets
std::vector<double> performAlgebraWithArrayFireGPU_batch(const std::vector<Eigen::VectorXd>& noise_batch, 
                                                        const std::vector<Eigen::VectorXd>& resvec_batch,
                                                        const std::vector<Eigen::VectorXd>& powercoeff_batch,
                                                        const std::vector<int>& totCoeff_batch,
                                                        const std::vector<double>& tdet_batch,
                                                        const std::vector<double>& freq_det_batch,
                                                        const std::vector<double>& timelike_batch,
                                                        const std::vector<double>& uniform_prior_batch);

// Debug function for comparing CPU and GPU implementations
// Note: Debug output has been removed for performance. This function is not called in production.
void compareEigenAndArrayFire(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                              double freq_det, double timelike, double uniform_prior);