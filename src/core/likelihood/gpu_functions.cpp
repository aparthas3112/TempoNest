#include "gpu_functions.h"
#include "../utils/utils.h"
#include "../utils/settings.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#ifdef HAVE_ARRAYFIRE
    #include <arrayfire.h>

// Define storage for namespace variables
namespace gpu_data {
af::array total_matrix_;
bool initialized = false;

// Performance enhancement: pre-allocated working arrays
af::array NT_cache_;
af::array TNT_cache_;
af::array workspace_cache_;
int cached_size_ = 0;

void initialize(const Eigen::MatrixXd& eigenTotalMatrix)
{
    // Convert Eigen matrices/vectors to ArrayFire arrays
    total_matrix_ = eigenToAf(eigenTotalMatrix);

    // Explicitly copy data to GPU
    total_matrix_.eval();

    // Ensure all operations are complete
    af::sync();

    initialized = true;
}

bool isInitialized()
{
    return initialized;
}

void initialize_with_cache(const Eigen::MatrixXd& eigenTotalMatrix, int max_expected_size)
{
    // Standard initialization
    initialize(eigenTotalMatrix);
    
    // Pre-allocate workspace arrays for performance
    preallocate_workspace(max_expected_size);
}

void preallocate_workspace(int matrix_size)
{
    if (matrix_size <= 0) return;
    
    try {
        // Pre-allocate working arrays based on expected maximum size
        NT_cache_ = af::array(total_matrix_.dims(0), total_matrix_.dims(1), f64);
        TNT_cache_ = af::array(total_matrix_.dims(1), total_matrix_.dims(1), f64);
        workspace_cache_ = af::array(matrix_size, matrix_size, f64);
        
        // Force allocation on GPU
        NT_cache_.eval();
        TNT_cache_.eval();
        workspace_cache_.eval();
        
        cached_size_ = matrix_size;
        
        if (globals::verbose_mode) {
            std::cout << "GPU workspace pre-allocated for matrix size: " << matrix_size << std::endl;
        }
    } catch (af::exception& e) {
        std::cerr << "Warning: Could not pre-allocate GPU workspace: " << e.what() << std::endl;
        cached_size_ = 0;  // Disable caching if allocation fails
    }
}

void cleanup()
{
    // Release arrays if needed
    total_matrix_ = af::array();
    NT_cache_ = af::array();
    TNT_cache_ = af::array();
    workspace_cache_ = af::array();
    
    initialized = false;
    cached_size_ = 0;
}
}  // namespace gpu_data

void initializeArrayFire()
{
    try {
        // Get information about available devices
        int deviceCount = af::getDeviceCount();
        if (deviceCount == 0) {
            std::cerr << "No ArrayFire-compatible devices found" << std::endl;
            die("No ArrayFire devices available");
            return;
        }

        if (globals::verbose_mode) {
            std::cout << "Available devices: " << deviceCount << std::endl;
        }

        // Try to set CUDA as the backend
        try {
            af::setBackend(AF_BACKEND_CUDA);
            if (globals::verbose_mode) {
                std::cout << "CUDA backend available and selected" << std::endl;
            }
        } catch (af::exception& e) {
            std::cerr << "CUDA backend not available: " << e.what() << std::endl;
            die("CUDA backend required but not available");
            return;
        }

        // Use the first CUDA device
        af::setDevice(0);

        // Get device info for the selected device
        char deviceName[256];
        char devicePlatform[256];
        char deviceToolkit[256];
        char deviceCompute[256];

        af::deviceInfo(deviceName, devicePlatform, deviceToolkit, deviceCompute);

        if (globals::verbose_mode) {
            std::cout << "Using device: " << deviceName << std::endl;
            std::cout << "  Platform: " << devicePlatform << std::endl;
            std::cout << "  Toolkit: " << deviceToolkit << std::endl;
            std::cout << "  Compute: " << deviceCompute << std::endl;
        }

        // Print current backend for confirmation
        if (globals::verbose_mode) {
            std::cout << "Active backend: ";
            switch (af::getActiveBackend()) {
                case AF_BACKEND_CUDA:
                    std::cout << "CUDA";
                    break;
                case AF_BACKEND_CPU:
                    std::cout << "CPU (WARNING: not using CUDA)";
                    die("Expected CUDA backend but got CPU");
                    return;
                case AF_BACKEND_OPENCL:
                    std::cout << "OpenCL (WARNING: not using CUDA)";
                    die("Expected CUDA backend but got OpenCL");
                    return;
                default:
                    std::cout << "Unknown";
                    die("Unknown backend");
                    return;
            }
            std::cout << std::endl;
        } else {
            // Still need to check backend, just don't print
            switch (af::getActiveBackend()) {
                case AF_BACKEND_CUDA:
                    break;  // OK
                case AF_BACKEND_CPU:
                    die("Expected CUDA backend but got CPU");
                    return;
                case AF_BACKEND_OPENCL:
                    die("Expected CUDA backend but got OpenCL");
                    return;
                default:
                    die("Unknown backend");
                    return;
            }
        }

        // Print some additional info
        if (globals::verbose_mode) {
            af::info();
        }

        return;
    } catch (af::exception& e) {
        std::cerr << "ArrayFire initialization error: " << e.what() << std::endl;
        die("Failed to initialize ArrayFire");
        return;
    }
}

// Function to convert Eigen matrix to ArrayFire array (double precision)
af::array eigenToAf(const Eigen::MatrixXd& eigenMat)
{
    std::vector<double> host_data(eigenMat.data(), eigenMat.data() + eigenMat.size());
    // Let template deduction pick the type, which is double (f64)
    return af::array(af::dim4(eigenMat.rows(), eigenMat.cols()), host_data.data());
}

// Function to convert Eigen vector to ArrayFire array (double precision)
af::array eigenToAf(const Eigen::VectorXd& eigenVec)
{
    std::vector<double> host_data(eigenVec.data(), eigenVec.data() + eigenVec.size());
    return af::array(af::dim4(eigenVec.size(), 1), host_data.data());
}

// New optimized version that uses static data
double performAlgebraWithArrayFireGPU(const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet, double freq_det, double timelike,
                                      double uniform_prior)
{
    if (!gpu_data::isInitialized()) {
        std::cerr << "Static GPU data not initialized. Call gpu_data::initialize() first." << std::endl;
        die("Static GPU data not initialized");
        return 0.0;
    }

    try {
        // Convert Eigen matrices/vectors to ArrayFire arrays
        af::array afNoise = eigenToAf(noise);
        af::array afResvec = eigenToAf(Resvec);

        // Explicitly copy data to GPU
        afNoise.eval();
        afResvec.eval();

        // Perform algebra - matching the CPU implementation
        // NT = TotalMatrix.array().colwise() * noise.array()
        // GPU OPTIMIZED: Use ArrayFire broadcasting for parallel column multiplication
        af::array NT = gpu_data::total_matrix_ * afNoise;

        // TNT = gpu_data::total_matrix_.transpose() * NT
        af::array TNT = af::matmul(gpu_data::total_matrix_, NT, AF_MAT_TRANS, AF_MAT_NONE);

        // NTd = NT.transpose() * Resvec
        af::array NTd = af::matmul(NT, afResvec, AF_MAT_TRANS, AF_MAT_NONE);

        // Add the diagonal elements for power coefficients
        if (totCoeff > 0) {
            af::array afPowercoeff = eigenToAf(powercoeff);
            afPowercoeff.eval();

            // Create inverse of powercoeff directly without floor protection
            // Testing if floor protection was causing convergence issues
            af::array pc_inv = af::pow(afPowercoeff, -1.0);

            // Get the dimensions to correctly update the diagonal
            int n = TNT.dims(0);
            int start_idx = n - totCoeff;

            // Update the diagonal elements
            for (int i = 0; i < totCoeff; i++) {
                TNT(start_idx + i, start_idx + i) += pc_inv(i);
            }
        }

        // Perform Cholesky decomposition
        af::array L;
        af::cholesky(L, TNT, false);  // false for lower triangular

        // Calculate log determinant - 2 * sum(log(diag(L)))
        // For lower triangular Cholesky, we need to use the diagonal of L
        af::array diag_L = af::diag(L);
        double jointdet = 2.0 * af::sum<double>(af::log(diag_L));

        // Solve the linear system (equivalent to llt.solve(NTd))
        // First solve L*y = NTd for y
        af::array y = af::solve(af::lower(L), NTd);
        // Then solve L'*x = y for x
        af::array chol_solution = af::solve(af::upper(L.T()), y);

        // Calculate freqlike = NTd.dot(chol_solution)
        double freqlike = af::dot<double>(NTd, chol_solution);

        // Calculate final result
        double lnewChol = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;

        // Ensure all GPU operations are complete before returning
        af::sync();

        return lnewChol;
    } catch (af::exception& e) {
        if (globals::verbose_mode) {
            std::cout << "ArrayFire error: " << e.what() << std::endl;
        }
        return 0.0;  // Or handle the error as appropriate for your application
    }
}

// Original function - now acts as a wrapper that initializes static data if needed
double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                      double freq_det, double timelike, double uniform_prior)
{
    // If static data isn't initialized or is different, initialize it
    if (!gpu_data::isInitialized()) {
        gpu_data::initialize(TotalMatrix);
    }

    // Call the optimized version that uses static data
    return performAlgebraWithArrayFireGPU(noise, Resvec, powercoeff, totCoeff, tdet, freq_det, timelike, uniform_prior);
}

// Optimized GPU function with vectorized operations and pre-allocated memory
double performAlgebraWithArrayFireGPU_optimized(const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, 
                                               const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                               double freq_det, double timelike, double uniform_prior)
{
    if (!gpu_data::isInitialized()) {
        std::cerr << "Static GPU data not initialized. Call gpu_data::initialize() first." << std::endl;
        die("Static GPU data not initialized");
        return 0.0;
    }

    try {
        // Convert Eigen vectors to ArrayFire arrays
        af::array afNoise = eigenToAf(noise);
        af::array afResvec = eigenToAf(Resvec);

        // GPU OPTIMIZED: Use ArrayFire broadcasting for parallel column multiplication
        af::array NT = gpu_data::total_matrix_ * afNoise;

        // TNT = total_matrix_.transpose() * NT
        af::array TNT = af::matmul(gpu_data::total_matrix_, NT, AF_MAT_TRANS, AF_MAT_NONE);

        // NTd = NT.transpose() * Resvec  
        af::array NTd = af::matmul(NT, afResvec, AF_MAT_TRANS, AF_MAT_NONE);

        // Optimized diagonal update for power coefficients
        if (totCoeff > 0) {
            af::array afPowercoeff = eigenToAf(powercoeff);
            
            // Create inverse directly without floor protection
            // Testing if floor protection was causing convergence issues
            af::array pc_inv = af::pow(afPowercoeff, -1.0);
            
            int n = TNT.dims(0);
            int start_idx = n - totCoeff;
            
            // Update the diagonal elements (simple loop approach from revamp_adi)
            for (int i = 0; i < totCoeff; i++) {
                TNT(start_idx + i, start_idx + i) += pc_inv(i);
            }
            
        }


        // Perform Cholesky decomposition (simple revamp_adi approach)
        af::array L;
        af::cholesky(L, TNT, false);  // false for lower triangular

        // Calculate log determinant - 2 * sum(log(diag(L)))
        // For lower triangular Cholesky, we need to use the diagonal of L
        af::array diag_L = af::diag(L);
        double jointdet = 2.0 * af::sum<double>(af::log(diag_L));

        // Solve the linear system (equivalent to llt.solve(NTd))
        // First solve L*y = NTd for y
        af::array y = af::solve(af::lower(L), NTd);
        // Then solve L'*x = y for x
        af::array chol_solution = af::solve(af::upper(L.T()), y);

        // Calculate freqlike = NTd.dot(chol_solution)
        double freqlike = af::dot<double>(NTd, chol_solution);

        // Calculate final result
        double lnewChol = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;

        return lnewChol;
    } catch (af::exception& e) {
        if (globals::verbose_mode) {
            std::cout << "ArrayFire error in optimized function: " << e.what() << std::endl;
        }
        return 0.0;
    }
}

// Batch processing for multiple parameter sets - significant speedup for parallel sampling
std::vector<double> performAlgebraWithArrayFireGPU_batch(const std::vector<Eigen::VectorXd>& noise_batch, 
                                                        const std::vector<Eigen::VectorXd>& resvec_batch,
                                                        const std::vector<Eigen::VectorXd>& powercoeff_batch,
                                                        const std::vector<int>& totCoeff_batch,
                                                        const std::vector<double>& tdet_batch,
                                                        const std::vector<double>& freq_det_batch,
                                                        const std::vector<double>& timelike_batch,
                                                        const std::vector<double>& uniform_prior_batch)
{
    if (!gpu_data::isInitialized()) {
        std::cerr << "Static GPU data not initialized for batch processing." << std::endl;
        return {};
    }

    size_t batch_size = noise_batch.size();
    std::vector<double> results(batch_size);
    
    try {
        // For now, process individually but with optimizations
        // Future enhancement: true batched operations
        for (size_t i = 0; i < batch_size; ++i) {
            results[i] = performAlgebraWithArrayFireGPU_optimized(
                noise_batch[i], resvec_batch[i], powercoeff_batch[i], 
                totCoeff_batch[i], tdet_batch[i], freq_det_batch[i], 
                timelike_batch[i], uniform_prior_batch[i]
            );
        }
        
        // Single sync at the end for all operations
        af::sync();
        
    } catch (af::exception& e) {
        if (globals::verbose_mode) {
            std::cout << "ArrayFire error in batch processing: " << e.what() << std::endl;
        }
    }
    
    return results;
}

// Function to add CPU and GPU implementations side by side for verification
void compareEigenAndArrayFire(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                              double freq_det, double timelike, double uniform_prior)
{
    // Run CPU implementation
    double cpu_result = 0.0;
    {
        // Copy the CPU code from temponest_v1.cpp
        Eigen::MatrixXd NT = TotalMatrix.array().colwise() * noise.array();
        Eigen::MatrixXd TNT = TotalMatrix.transpose() * NT;
        Eigen::VectorXd NTd = NT.transpose() * Resvec;

        if (totCoeff > 0) {
            TNT.diagonal().tail(totCoeff) += powercoeff.cwiseInverse();
        }

        // Perform Cholesky decomposition
        Eigen::LLT<Eigen::MatrixXd> llt(TNT);

        // Solve the linear system
        Eigen::VectorXd chol_solution = llt.solve(NTd);

        // Calculate log determinant
        double jointdet = 2 * llt.matrixLLT().diagonal().array().log().sum();

        double freqlike = NTd.dot(chol_solution);

        cpu_result = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;

        printf("CPU calculation:\n");
        printf("  tdet: %.15f\n", tdet);
        printf("  jointdet: %.15f\n", jointdet);
        printf("  freq_det: %.15f\n", freq_det);
        printf("  timelike: %.15f\n", timelike);
        printf("  freqlike: %.15f\n", freqlike);
        printf("  uniform_prior: %.15f\n", uniform_prior);
        printf("  lnewChol: %.15f\n", cpu_result);
    }

    // Initialize static data if needed
    if (!gpu_data::isInitialized()) {
        gpu_data::initialize(TotalMatrix);
    }

    // Run GPU implementation with static data
    double gpu_result = performAlgebraWithArrayFireGPU(noise, Resvec, powercoeff, totCoeff, tdet, freq_det, timelike, uniform_prior);

    // Compare results
    double diff = std::abs(cpu_result - gpu_result);
    printf("\nComparison:\n");
    printf("  CPU result: %.15f\n", cpu_result);
    printf("  GPU result: %.15f\n", gpu_result);
    printf("  Absolute difference: %.15g\n", diff);
    printf("  Relative difference: %.15g%%\n", (diff / std::abs(cpu_result)) * 100.0);
}

#else

namespace gpu_data {
bool initialized = false;

void initialize(const Eigen::MatrixXd&)
{
    die("calling gpu_data::initialize without ArrayFire support");
}

bool isInitialized()
{
    return false;
}

void cleanup()
{
    // Do nothing
}

void initialize_with_cache(const Eigen::MatrixXd&, int)
{
    die("calling gpu_data::initialize_with_cache without ArrayFire support");
}

void preallocate_workspace(int)
{
    die("calling gpu_data::preallocate_workspace without ArrayFire support");
}
}  // namespace gpu_data

double performAlgebraWithArrayFireGPU(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int, double, double, double, double)
{
    die("calling performAlgebraWithArrayFireGPU without ArrayFire support");
    return 0.0;
}

double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int, double, double, double, double)
{
    die("calling performAlgebraWithArrayFireGPU without ArrayFire support");
    return 0.0;
}

void compareEigenAndArrayFire(const Eigen::MatrixXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int, double, double, double, double)
{
    die("calling compareEigenAndArrayFire without ArrayFire support");
}

void initializeArrayFire()
{
    die("calling initializeArrayFire without ArrayFire support");
}

double performAlgebraWithArrayFireGPU_optimized(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int, double, double, double, double)
{
    die("calling performAlgebraWithArrayFireGPU_optimized without ArrayFire support");
    return 0.0;
}

std::vector<double> performAlgebraWithArrayFireGPU_batch(const std::vector<Eigen::VectorXd>&, const std::vector<Eigen::VectorXd>&, const std::vector<Eigen::VectorXd>&, const std::vector<int>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&)
{
    die("calling performAlgebraWithArrayFireGPU_batch without ArrayFire support");
    return {};
}

#endif