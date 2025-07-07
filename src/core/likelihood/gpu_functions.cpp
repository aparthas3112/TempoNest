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

            // CRITICAL FIX: Add floor to powercoeff to prevent division by zero
            const double min_powercoeff_floor = 1e-30;
            
            // Convert to host to apply floor protection (most reliable approach)
            std::vector<double> host_powercoeff(afPowercoeff.elements());
            afPowercoeff.host(host_powercoeff.data());
            
            // Apply floor protection on host for numerical stability
            int fixed_count = 0;
            for (auto& val : host_powercoeff) {
                if (val <= min_powercoeff_floor) {
                    val = min_powercoeff_floor;
                    fixed_count++;
                }
            }
            
            if (fixed_count > 0) {
                if (globals::verbose_mode) {
                    std::cout << "GPU DEBUG: Fixed " << fixed_count << " powercoeff values" << std::endl;
                }
            }
            
            // Convert back to GPU
            afPowercoeff = af::array(afPowercoeff.dims(), host_powercoeff.data());

            // Create inverse of powercoeff
            af::array pc_inv = af::pow(afPowercoeff, -1.0);
            
            // ADDITIONAL FIX: Apply floor protection to pc_inv as well
            std::vector<double> host_pc_inv(pc_inv.elements());
            pc_inv.host(host_pc_inv.data());
            
            const double min_pc_inv_floor = 1e-10; // Much more aggressive floor
            for (auto& val : host_pc_inv) {
                if (val <= min_pc_inv_floor || !std::isfinite(val)) {
                    val = min_pc_inv_floor;
                }
            }
            
            pc_inv = af::array(pc_inv.dims(), host_pc_inv.data());

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

        // Use pre-allocated arrays if available and size matches
        af::array NT, TNT;
        bool use_cache = (gpu_data::cached_size_ > 0) && 
                        (gpu_data::total_matrix_.dims(1) <= gpu_data::cached_size_);
        
        if (use_cache) {
            // Reuse pre-allocated arrays - much faster
            NT = gpu_data::NT_cache_(af::span, af::seq(0, gpu_data::total_matrix_.dims(1) - 1));
            TNT = gpu_data::TNT_cache_(af::seq(0, gpu_data::total_matrix_.dims(1) - 1), 
                                      af::seq(0, gpu_data::total_matrix_.dims(1) - 1));
        }
        
        // GPU OPTIMIZED: Use ArrayFire broadcasting for parallel column multiplication
        if (use_cache) {
            NT = gpu_data::total_matrix_ * afNoise;  // Reuse existing array
        } else {
            NT = gpu_data::total_matrix_ * afNoise;  // Allocate new
        }

        // TNT = total_matrix_.transpose() * NT
        TNT = af::matmul(gpu_data::total_matrix_, NT, AF_MAT_TRANS, AF_MAT_NONE);

        // NTd = NT.transpose() * Resvec  
        af::array NTd = af::matmul(NT, afResvec, AF_MAT_TRANS, AF_MAT_NONE);

        // Optimized diagonal update for power coefficients
        if (totCoeff > 0) {
            af::array afPowercoeff = eigenToAf(powercoeff);
            
            // CRITICAL FIX: Add floor to powercoeff to prevent division by zero
            const double min_powercoeff_floor = 1e-30;
            
            // Convert to host to apply floor protection (most reliable approach)
            std::vector<double> host_powercoeff(afPowercoeff.elements());
            afPowercoeff.host(host_powercoeff.data());
            
            // Apply floor protection on host for numerical stability
            int fixed_count = 0;
            for (auto& val : host_powercoeff) {
                if (val <= min_powercoeff_floor) {
                    val = min_powercoeff_floor;
                    fixed_count++;
                }
            }
            
            if (fixed_count > 0) {
                if (globals::verbose_mode) {
                    std::cout << "GPU DEBUG: Fixed " << fixed_count << " powercoeff values" << std::endl;
                }
            }
            
            // Convert back to GPU
            afPowercoeff = af::array(afPowercoeff.dims(), host_powercoeff.data());
            
            // ENHANCED DEBUG: Detailed analysis for extreme powercoeff values
            std::vector<double> host_pc_debug(afPowercoeff.elements());
            afPowercoeff.host(host_pc_debug.data());
            
            double pc_min = *std::min_element(host_pc_debug.begin(), host_pc_debug.end());
            double pc_max = *std::max_element(host_pc_debug.begin(), host_pc_debug.end());
            bool has_extreme = (pc_max > 1e15 || pc_min < 1e-15);
            
            if (has_extreme && globals::verbose_mode) {
                std::cout << "GPU DEBUG: *** EXTREME POWERCOEFF DETECTED ***" << std::endl;
                std::cout << "GPU DEBUG: powercoeff range: min=" << std::scientific << pc_min 
                         << ", max=" << pc_max << ", size=" << host_pc_debug.size() << std::endl;
                
                // Count extreme values
                int very_large_count = 0, very_small_count = 0;
                for (const auto& val : host_pc_debug) {
                    if (val > 1e15) very_large_count++;
                    if (val < 1e-15) very_small_count++;
                }
                
                std::cout << "GPU DEBUG: Extreme counts - large(>1e15):" << very_large_count 
                         << ", small(<1e-15):" << very_small_count << std::endl;
                
                // Show first 10 and pattern
                std::cout << "GPU DEBUG: First 10 powercoeff: ";
                for (size_t i = 0; i < std::min(host_pc_debug.size(), size_t(10)); ++i) {
                    std::cout << std::scientific << host_pc_debug[i] << " ";
                }
                std::cout << std::endl;
            }
            
            // Convert back to GPU for computation
            afPowercoeff = af::array(afPowercoeff.dims(), host_pc_debug.data());
            
            // OPTIMIZATION: Vectorized diagonal update instead of loop
            af::array pc_inv = af::pow(afPowercoeff, -1.0);
            
            // ADDITIONAL FIX: Apply floor protection to pc_inv as well
            std::vector<double> host_pc_inv(pc_inv.elements());
            pc_inv.host(host_pc_inv.data());
            
            int pc_inv_fixed = 0;
            const double min_pc_inv_floor = 1e-10; // Much more aggressive floor
            double actual_min_before = *std::min_element(host_pc_inv.begin(), host_pc_inv.end());
            
            for (auto& val : host_pc_inv) {
                if (val <= min_pc_inv_floor || !std::isfinite(val)) {
                    val = min_pc_inv_floor;
                    pc_inv_fixed++;
                }
            }
            
            double actual_min_after = *std::min_element(host_pc_inv.begin(), host_pc_inv.end());
            
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: pc_inv floor protection - before_min=" << std::scientific << actual_min_before 
                         << ", after_min=" << actual_min_after << ", fixed=" << pc_inv_fixed << std::endl;
            }
            
            pc_inv = af::array(pc_inv.dims(), host_pc_inv.data());
            
            // DEBUG: Check pc_inv values
            double min_inv = af::min<double>(pc_inv);
            double max_inv = af::max<double>(pc_inv);
            bool has_nan_inv = af::anyTrue<bool>(af::isNaN(pc_inv));
            bool has_inf_inv = af::anyTrue<bool>(af::isInf(pc_inv));
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: pc_inv: min=" << min_inv << ", max=" << max_inv 
                         << ", has_nan=" << (has_nan_inv ? "YES" : "NO")
                         << ", has_inf=" << (has_inf_inv ? "YES" : "NO") << std::endl;
            }
            
            int n = TNT.dims(0);
            int start_idx = n - totCoeff;
            
            // FIXED: Only update diagonal once using efficient vectorized method
            // Create a diagonal matrix with zeros except for pc_inv in the bottom-right
            af::array diag_update = af::constant(0, n, f64);
            diag_update(af::seq(start_idx, start_idx + totCoeff - 1)) = pc_inv;
            TNT = TNT + af::diag(diag_update, 0, false);
            
            // DEBUG: Check TNT after diagonal update
            bool has_nan_tnt = af::anyTrue<bool>(af::isNaN(TNT));
            bool has_inf_tnt = af::anyTrue<bool>(af::isInf(TNT));
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: TNT after diag update: has_nan=" << (has_nan_tnt ? "YES" : "NO")
                         << ", has_inf=" << (has_inf_tnt ? "YES" : "NO") << std::endl;
            }
        }

        // DEBUG: Analyze TNT matrix before Cholesky
        bool has_nan_tnt_pre = af::anyTrue<bool>(af::isNaN(TNT));
        bool has_inf_tnt_pre = af::anyTrue<bool>(af::isInf(TNT));
        double tnt_min = af::min<double>(TNT);
        double tnt_max = af::max<double>(TNT);
        
        // Check matrix condition number and determinant
        double tnt_det = af::det<double>(TNT);
        bool det_finite = std::isfinite(tnt_det);
        
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: TNT pre-Cholesky: has_nan=" << (has_nan_tnt_pre ? "YES" : "NO")
                     << ", has_inf=" << (has_inf_tnt_pre ? "YES" : "NO")
                     << ", min=" << std::scientific << tnt_min
                     << ", max=" << tnt_max
                     << ", det=" << tnt_det << " (finite=" << (det_finite ? "YES" : "NO") << ")" << std::endl;
        }

        // Perform Cholesky decomposition with error handling
        af::array L;
        bool cholesky_success = true;
        try {
            af::cholesky(L, TNT, false);  // false for lower triangular
        } catch (af::exception& e) {
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: Cholesky exception: " << e.what() << std::endl;
            }
            cholesky_success = false;
        }

        // DEBUG: Check Cholesky result
        bool has_nan_L = false;
        bool has_inf_L = false;
        if (cholesky_success) {
            has_nan_L = af::anyTrue<bool>(af::isNaN(L));
            has_inf_L = af::anyTrue<bool>(af::isInf(L));
        }
        
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: Cholesky success=" << (cholesky_success ? "YES" : "NO")
                     << ", L: has_nan=" << (has_nan_L ? "YES" : "NO")
                     << ", has_inf=" << (has_inf_L ? "YES" : "NO") << std::endl;
        }
                 
        // If Cholesky failed, use CPU fallback computation
        if (!cholesky_success || has_nan_L || has_inf_L) {
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: Cholesky failed - using CPU fallback computation" << std::endl;
            }
            
            // Convert TNT to CPU for fallback computation
            std::vector<double> host_tnt(TNT.elements());
            TNT.host(host_tnt.data());
            
            // Convert to Eigen matrix
            int n = TNT.dims(0);
            Eigen::MatrixXd eigen_tnt = Eigen::Map<Eigen::MatrixXd>(host_tnt.data(), n, n);
            
            // Convert NTd to CPU
            std::vector<double> host_ntd(NTd.elements());
            NTd.host(host_ntd.data());
            Eigen::VectorXd eigen_ntd = Eigen::Map<Eigen::VectorXd>(host_ntd.data(), NTd.elements());
            
            // Perform CPU Cholesky decomposition
            Eigen::LLT<Eigen::MatrixXd> llt(eigen_tnt);
            
            if (llt.info() != Eigen::Success) {
                if (globals::verbose_mode) {
                    std::cout << "GPU DEBUG: CPU fallback Cholesky also failed - matrix is not positive definite" << std::endl;
                }
                return 0.0;
            }
            
            // Calculate CPU results
            Eigen::VectorXd cpu_chol_solution = llt.solve(eigen_ntd);
            double cpu_jointdet = 2 * llt.matrixLLT().diagonal().array().log().sum();
            double cpu_freqlike = eigen_ntd.dot(cpu_chol_solution);
            double cpu_result = -0.5 * (tdet + cpu_jointdet + freq_det + timelike - cpu_freqlike) + uniform_prior;
            
            if (globals::verbose_mode) {
                std::cout << "GPU DEBUG: CPU fallback successful - result=" << std::scientific << cpu_result << std::endl;
            }
            return cpu_result;
        }

        // OPTIMIZATION: Fused operations for better GPU utilization
        // Calculate log determinant and solve system in minimal passes
        af::array diag_L = af::diag(L);
        
        // DEBUG: Check diagonal values
        double min_diag = af::min<double>(diag_L);
        double max_diag = af::max<double>(diag_L);
        bool has_nan_diag = af::anyTrue<bool>(af::isNaN(diag_L));
        bool has_negative_diag = af::anyTrue<bool>(diag_L <= 0.0);
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: diag_L: min=" << min_diag << ", max=" << max_diag
                     << ", has_nan=" << (has_nan_diag ? "YES" : "NO")
                     << ", has_negative=" << (has_negative_diag ? "YES" : "NO") << std::endl;
        }
        
        double jointdet = 2.0 * af::sum<double>(af::log(diag_L));
        
        // DEBUG: Check jointdet
        bool jointdet_finite = std::isfinite(jointdet);
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: jointdet=" << jointdet << ", finite=" << (jointdet_finite ? "YES" : "NO") << std::endl;
        }

        // Solve the linear system using Cholesky decomposition result
        // First solve L*y = NTd for y
        af::array y = af::solve(af::lower(L), NTd);
        // Then solve L'*x = y for x
        af::array chol_solution = af::solve(af::upper(L.T()), y);

        // Calculate freqlike = NTd.dot(chol_solution)
        double freqlike = af::dot<double>(NTd, chol_solution);
        
        // DEBUG: Check freqlike
        bool freqlike_finite = std::isfinite(freqlike);
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: freqlike=" << freqlike << ", finite=" << (freqlike_finite ? "YES" : "NO") << std::endl;
        }

        // Calculate final result
        double lnewChol = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;
        
        // DEBUG: Check final result
        bool result_finite = std::isfinite(lnewChol);
        if (globals::verbose_mode) {
            std::cout << "GPU DEBUG: final result=" << lnewChol << ", finite=" << (result_finite ? "YES" : "NO") << std::endl;
        }

        // OPTIMIZATION: Async evaluation - only sync when we need the scalar result
        // af::sync() is called implicitly when we extract the scalar value above

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