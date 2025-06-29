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

// Initialize the static matrices
void initialize(const Eigen::MatrixXd& eigenTotalMatrix);

// Check if initialized
bool isInitialized();

// Clean up resources if needed
void cleanup();
}  // namespace gpu_data

#else
// Stub namespace when ArrayFire is not available
namespace gpu_data {
extern bool initialized;

void initialize(const Eigen::MatrixXd&);
bool isInitialized();
void cleanup();
}  // namespace gpu_data
#endif

void initializeArrayFire();
double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                      double freq_det, double timelike, double uniform_prior);
void compareEigenAndArrayFire(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                              double freq_det, double timelike, double uniform_prior);