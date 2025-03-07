#include <utils.h>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#ifdef HAVE_ARRAYFIRE
    #include <arrayfire.h>

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

        std::cout << "Available devices: " << deviceCount << std::endl;

        // Try to set CUDA as the backend
        try {
            af::setBackend(AF_BACKEND_CUDA);
            std::cout << "CUDA backend available and selected" << std::endl;
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

        std::cout << "Using device: " << deviceName << std::endl;
        std::cout << "  Platform: " << devicePlatform << std::endl;
        std::cout << "  Toolkit: " << deviceToolkit << std::endl;
        std::cout << "  Compute: " << deviceCompute << std::endl;

        // Print current backend for confirmation
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

        // Print some additional info
        af::info();

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

double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                      double freq_det, double timelike, double uniform_prior)
{

    try {
        // Convert Eigen matrices/vectors to ArrayFire arrays
        af::array afTotalMatrix = eigenToAf(TotalMatrix);
        af::array afNoise = eigenToAf(noise);
        af::array afResvec = eigenToAf(Resvec);

        // Explicitly copy data to GPU
        afTotalMatrix.eval();
        afNoise.eval();
        afResvec.eval();

        // Verify array precision (for debugging)
        // std::cout << "TotalMatrix dtype: " << afTotalMatrix.type() << std::endl;
        // std::cout << "noise dtype: " << afNoise.type() << std::endl;
        // std::cout << "Resvec dtype: " << afResvec.type() << std::endl;
        // af::dtype::f64 should be displayed as 2

        // Perform algebra - matching the CPU implementation
        // NT = TotalMatrix.array().colwise() * noise.array()
        // We need to multiply each column of TotalMatrix by the corresponding element in noise
        af::array NT = af::constant(0, afTotalMatrix.dims(), f64);
        for (int i = 0; i < TotalMatrix.cols(); i++) {
            NT(af::span, i) = afTotalMatrix(af::span, i) * afNoise;
        }

        // TNT = TotalMatrix.transpose() * NT
        af::array TNT = af::matmul(afTotalMatrix, NT, AF_MAT_TRANS, AF_MAT_NONE);

        // NTd = NT.transpose() * Resvec
        af::array NTd = af::matmul(NT, afResvec, AF_MAT_TRANS, AF_MAT_NONE);

        // Add the diagonal elements for power coefficients
        if (totCoeff > 0) {
            af::array afPowercoeff = eigenToAf(powercoeff);
            afPowercoeff.eval();

            // Create inverse of powercoeff
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

        // Debug output
        /*
        printf("GPU calculation:\n");
        printf("  tdet: %.15f\n", tdet);
        printf("  jointdet: %.15f\n", jointdet);
        printf("  freq_det: %.15f\n", freq_det);
        printf("  timelike: %.15f\n", timelike);
        printf("  freqlike: %.15f\n", freqlike);
        printf("  uniform_prior: %.15f\n", uniform_prior);
        printf("  lnewChol: %.15f\n", lnewChol);
        */
        // Ensure all GPU operations are complete before returning
        af::sync();

        return lnewChol;
    } catch (af::exception& e) {
        std::cout << "ArrayFire error: " << e.what() << std::endl;
        return 0.0;  // Or handle the error as appropriate for your application
    }
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

    // Run GPU implementation
    double gpu_result = performAlgebraWithArrayFireGPU(TotalMatrix, noise, Resvec, powercoeff, totCoeff, tdet, freq_det, timelike, uniform_prior);

    // Compare results
    double diff = std::abs(cpu_result - gpu_result);
    printf("\nComparison:\n");
    printf("  CPU result: %.15f\n", cpu_result);
    printf("  GPU result: %.15f\n", gpu_result);
    printf("  Absolute difference: %.15g\n", diff);
    printf("  Relative difference: %.15g%%\n", (diff / std::abs(cpu_result)) * 100.0);
}

#else

double performAlgebraWithArrayFireGPU(const Eigen::MatrixXd& TotalMatrix, const Eigen::VectorXd& noise, const Eigen::VectorXd& Resvec, const Eigen::VectorXd& powercoeff, int totCoeff, double tdet,
                                      double freq_det, double timelike, double uniform_prior)
{
    die("calling performAlgebraWithArrayFireGPU without ArrayFire support");
    return 0.0;
}

void initializeArrayFire()
{
    die("calling initializeArrayFire without ArrayFire support");
}

#endif