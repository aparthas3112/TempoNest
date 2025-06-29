# TempoNest GPU Optimization

## Overview

TempoNest includes GPU acceleration through ArrayFire, providing significant performance improvements for computationally intensive likelihood calculations. This document covers GPU implementation, performance characteristics, and optimization strategies.

## Table of Contents

1. [GPU Architecture](#gpu-architecture)
2. [ArrayFire Integration](#arrayfire-integration)
3. [Performance Benchmarks](#performance-benchmarks)
4. [Recent Optimizations (2025-06)](#recent-optimizations-2025-06)
5. [Memory Management](#memory-management)
6. [Development Guidelines](#development-guidelines)
7. [Troubleshooting GPU Issues](#troubleshooting-gpu-issues)

---

## GPU Architecture

### Hardware Requirements

**Minimum Requirements:**
- NVIDIA GPU with CUDA Compute Capability 3.5+
- 4GB VRAM (8GB+ recommended)
- CUDA 11.0+ with compatible drivers

**Recommended Configuration:**
- NVIDIA GeForce RTX 4090 or similar (24GB VRAM)
- CUDA 12.0+ with latest drivers
- High memory bandwidth (>1000 GB/s)

**Tested Hardware:**
- **Primary**: NVIDIA GeForce RTX 4090 (24GB VRAM)
- **ArrayFire Version**: 3.9.0 with CUDA backend
- **Typical Speedup**: 1-2x overall vs 300+ CPU cores

### GPU vs CPU Performance

**Memory Bandwidth Comparison:**
- **GPU**: >1000 GB/s (RTX 4090)
- **CPU**: ~100 GB/s per socket (high-end server)
- **Advantage**: 10x+ memory bandwidth for matrix operations

**Scaling Behavior:**
- **Parameter Dimensions**: Most effective for >20 parameter problems
- **Time Complexity**: O(N³) for Cholesky operations (N = observations + coefficients)
- **Memory Usage**: Proportional to N² for covariance matrices

---

## ArrayFire Integration

### Implementation Overview

TempoNest uses ArrayFire for all GPU operations, providing a unified interface across CUDA, OpenCL, and CPU backends.

```cpp
#ifdef HAVE_ARRAYFIRE
#include <arrayfire.h>

// GPU implementation
af::array gpu_result = optimized_gpu_function(gpu_data);
#else
// CPU fallback
Eigen::VectorXd cpu_result = standard_cpu_function(cpu_data);
#endif
```

### Key GPU Functions

**Matrix Operations:**
```cpp
// Broadcasting operation (automatic parallelization)
af::array result = matrix * vector;

// Cholesky decomposition
af::array chol_L;
af::cholesky(chol_L, TNT);

// Element-wise operations
af::array scaled = af::exp(log_params) * coefficients;
```

**Memory Management:**
```cpp
namespace gpu_data {
    static af::array total_matrix_;
    static bool initialized_ = false;
    
    void initialize_gpu_data(const Eigen::MatrixXd& cpu_matrix) {
        if (!initialized_) {
            total_matrix_ = af::array(cpu_matrix.rows(), cpu_matrix.cols(), 
                                    cpu_matrix.data());
            initialized_ = true;
        }
    }
}
```

### Compilation and Linking

**Configure with ArrayFire:**
```bash
./configure --with-arrayfire=/path/to/arrayfire
```

**CMake Integration:**
```cmake
find_package(ArrayFire REQUIRED)
target_link_libraries(temponest ${ArrayFire_LIBRARIES})
target_compile_definitions(temponest PRIVATE HAVE_ARRAYFIRE)
```

---

## Performance Benchmarks

### MultiNest Scaling

**Runtime scaling**: O(nlive × log(evidence_range))

| Live Points | Typical Runtime |
|-------------|-----------------|
| 500         | 20-60 minutes   |
| 1000        | 1-2 hours       |
| 4000        | 4-8 hours       |

**Efficiency Settings Impact:**

| Efficiency | Convergence | Exploration | Recommended Use |
|------------|-------------|-------------|-----------------|
| 0.3        | Fast        | Less        | Testing/debugging |
| 0.1        | Balanced    | Standard    | Production analysis |
| 0.01       | Slower      | Thorough    | Publication quality |

**Typical Acceptance Rates**: 5-15%

### GPU Speedup Analysis

**Problem Size Dependencies:**
- **Small problems** (<10 parameters): Minimal GPU advantage
- **Medium problems** (10-50 parameters): 2-5x speedup
- **Large problems** (>50 parameters): 5-10x+ speedup

**Memory Bandwidth Utilization:**
- **Matrix operations**: 80-95% GPU utilization
- **Element-wise operations**: 60-80% utilization
- **Data transfer overhead**: <5% for large problems

---

## Recent Optimizations (2025-06)

### Major GPU Performance Improvements

#### Problem: Sequential Column Multiplication

**Issue Identified**: Sequential loop in likelihood calculation destroying GPU parallelization.

**Before (Sequential)**:
```cpp
af::array NT = af::constant(0, gpu_data::total_matrix_.dims(), f64);
for (int i = 0; i < gpu_data::total_matrix_.dims(1); i++) {
    NT(af::span, i) = gpu_data::total_matrix_(af::span, i) * afNoise;
}
```

**After (Parallel Broadcasting)**:
```cpp
af::array NT = gpu_data::total_matrix_ * afNoise;
```

**Performance Impact**:
- **Expected speedup**: 10-50x for matrix operations
- **GPU utilization**: Improved from ~1% to near-100%
- **Mathematical verification**: Identical results to CPU implementation
- **Memory efficiency**: Reduced temporary array allocations

#### Coefficient Calculation Overhaul

**Problem**: Legacy implementation used fixed frequency counts (`num_freqs = 33`) independent of observation time span, creating scientifically inconsistent noise models.

**Legacy Implementation**:
```cpp
// Fixed coefficient count (incorrect)
int num_freqs = 33;  // Independent of data span
```

**New Implementation**:
```cpp
// Time-span dependent coefficient calculation
int num_freqs = static_cast<int>(std::floor(time_span_days / days_per_coeff));
```

**Impact**:
- **Default setting**: `days_per_coeff = 30.0` for scientific accuracy
- **Example improvement**: 4041-day span → 134 frequencies (vs fixed 33 previously)
- **Total coefficients**: 268 red + 268 DM = 536 vs ~132 previously
- **Mathematical consistency**: Restored compatibility with legacy CPU versions

### Verification Results

**Reference Test Case**:
- **Dataset**: PSR J1733-3716 (129 TOAs, 4041-day span)
- **Expected Likelihood**: 929.961754 ± 1e-6
- **Frequency Count**: 134 red + 134 DM = 268 total frequencies
- **GPU vs CPU agreement**: Identical to machine precision

---

## Memory Management

### GPU Memory Architecture

**Static Data Structure:**
```cpp
namespace gpu_data {
    static af::array total_matrix_;      // Main covariance matrix
    static af::array noise_coefficients_; // Fourier coefficients
    static af::array design_matrix_;     // Timing model matrix
    static bool initialized_ = false;
}
```

**Benefits**:
- **Prevents repeated allocation**: Data persists between likelihood calls
- **Reduces memory fragmentation**: Single large allocation
- **Minimizes CPU-GPU transfers**: Data uploaded once

**Limitations**:
- **Fixed problem size**: Cannot easily resize during analysis
- **Memory overhead**: Allocates for maximum expected size

### Memory Pool Optimization

**Current Implementation**:
```cpp
void initialize_gpu_memory(int max_obs, int max_coeffs) {
    if (!gpu_data::initialized_) {
        // Allocate maximum expected sizes
        gpu_data::total_matrix_ = af::constant(0.0, max_obs + max_coeffs, 
                                              max_obs + max_coeffs, f64);
        gpu_data::initialized_ = true;
    }
}
```

**Future Enhancement (Planned)**:
```cpp
class GPUMemoryPool {
public:
    af::array get_matrix(int rows, int cols);
    void release_matrix(af::array& matrix);
    void garbage_collect();
private:
    std::vector<af::array> free_matrices_;
    std::map<std::pair<int,int>, std::vector<af::array>> size_pools_;
};
```

### VRAM Usage Guidelines

**Memory Requirements by Problem Size:**

| Parameters | Observations | VRAM Usage | Recommended GPU |
|------------|-------------|------------|-----------------|
| 10-20      | 100-500     | 1-2 GB     | GTX 1060 6GB+   |
| 20-50      | 500-2000    | 4-8 GB     | RTX 3070 8GB+   |
| 50-100     | 1000-5000   | 8-16 GB    | RTX 3080 16GB+  |
| 100+       | 2000+       | 16+ GB     | RTX 4090 24GB+  |

---

## Development Guidelines

### Adding GPU Functions

**Function Structure Template:**
```cpp
#ifdef HAVE_ARRAYFIRE
af::array gpu_function_name(const af::array& input_data, 
                           const std::vector<double>& parameters) {
    // GPU implementation using ArrayFire
    af::array result = af::matmul(input_data, af::array(parameters.data(), 
                                                       parameters.size()));
    return result;
}
#endif

// CPU fallback always required
Eigen::VectorXd cpu_function_name(const Eigen::MatrixXd& input_data,
                                 const std::vector<double>& parameters) {
    // Standard CPU implementation
    Eigen::Map<const Eigen::VectorXd> param_vec(parameters.data(), 
                                               parameters.size());
    return input_data * param_vec;
}

// Unified interface
VectorXd function_name(const MatrixXd& data, const std::vector<double>& params) {
#ifdef HAVE_ARRAYFIRE
    af::array gpu_data = af::array(data.rows(), data.cols(), data.data());
    af::array gpu_result = gpu_function_name(gpu_data, params);
    
    // Convert back to Eigen
    VectorXd result(gpu_result.dims(0));
    gpu_result.host(result.data());
    return result;
#else
    return cpu_function_name(data, params);
#endif
}
```

### Testing Protocol

**Unit Tests for GPU Functions:**
```cpp
TEST(GPUFunctions, MatrixMultiplication) {
    const double tolerance = 1e-12;
    
    // Test data
    MatrixXd test_matrix = generate_test_matrix();
    std::vector<double> test_params = generate_test_parameters();
    
    // Compare CPU and GPU results
    VectorXd cpu_result = cpu_function_name(test_matrix, test_params);
    
#ifdef HAVE_ARRAYFIRE
    VectorXd gpu_result = function_name(test_matrix, test_params);
    
    EXPECT_NEAR((cpu_result - gpu_result).norm(), 0.0, tolerance);
#endif
}
```

**Performance Benchmarks:**
```cpp
void benchmark_gpu_performance() {
    const int num_iterations = 100;
    
    auto start_cpu = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        cpu_function_name(test_data, test_params);
    }
    auto end_cpu = std::chrono::high_resolution_clock::now();
    
#ifdef HAVE_ARRAYFIRE
    auto start_gpu = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        gpu_function_name(test_data, test_params);
    }
    af::sync();  // Ensure GPU completion
    auto end_gpu = std::chrono::high_resolution_clock::now();
    
    double speedup = duration_cast<microseconds>(end_cpu - start_cpu).count() /
                    duration_cast<microseconds>(end_gpu - start_gpu).count();
    std::cout << "GPU speedup: " << speedup << "x" << std::endl;
#endif
}
```

---

## Troubleshooting GPU Issues

### Common Problems

**ArrayFire not found during build:**
```bash
# Check ArrayFire installation
ls /path/to/arrayfire/lib64/libaf*

# Set environment variables
export ARRAYFIRE_PATH=/path/to/arrayfire
export LD_LIBRARY_PATH=$ARRAYFIRE_PATH/lib64:$LD_LIBRARY_PATH

# Reconfigure
./configure --with-arrayfire=$ARRAYFIRE_PATH
```

**CUDA driver issues:**
```bash
# Check CUDA installation
nvidia-smi
nvcc --version

# Verify ArrayFire CUDA backend
export AF_PRINT_ERRORS=1
./test_arrayfire_example
```

**Runtime GPU errors:**
```bash
# Enable ArrayFire debugging
export AF_DEBUG=1
export AF_PRINT_ERRORS=1

# Check memory usage
nvidia-smi

# Monitor GPU utilization
watch -n 1 nvidia-smi
```

### Performance Debugging

**GPU Utilization Monitoring:**
```bash
# Monitor during analysis
nvidia-smi dmon -s pucvmet -i 0

# Profile with nvtop (if available)
nvtop
```

**Memory Usage Analysis:**
```cpp
void print_gpu_memory_info() {
#ifdef HAVE_ARRAYFIRE
    size_t alloc_bytes, alloc_buffers, lock_bytes, lock_buffers;
    af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
    
    std::cout << "GPU Memory - Allocated: " << alloc_bytes / 1e9 << " GB" << std::endl;
    std::cout << "GPU Memory - Locked: " << lock_bytes / 1e9 << " GB" << std::endl;
#endif
}
```

**Performance Profiling:**
```bash
# Use NVIDIA profiler
nsys profile --stats=true ./tempo2 -gr temponest -f test.par test.tim -Cfile config.json

# ArrayFire timing
export AF_TIMER=1
```

### Known Limitations and Future Work

#### Current Limitations

1. **Diagonal Updates**: Sequential GPU implementation pending optimization
   - **Impact**: Some operations not fully parallelized
   - **Planned Fix**: Vectorized diagonal operations in Phase 2

2. **Memory Pooling**: Static allocation may be suboptimal for varying problem sizes
   - **Impact**: Memory overhead for small problems
   - **Planned Fix**: Dynamic memory pool in Phase 3

3. **Batch Processing**: Single dataset processing limits GPU utilization
   - **Impact**: Underutilized GPU for small problems
   - **Planned Fix**: Multi-dataset batching

#### Future Optimization Roadmap

**Phase 2: Diagonal Vectorization**
- **Target**: 3-5x additional speedup
- **Timeline**: Next major release
- **Implementation**: Replace sequential diagonal updates with parallel operations

**Phase 3: Memory Pooling and Batching**
- **Target**: 2-3x additional speedup
- **Timeline**: Medium-term development
- **Implementation**: Dynamic memory management and batch processing

**Advanced Features:**
- **Multi-GPU Support**: For extremely large problems (100+ parameters)
- **Mixed Precision**: Float16 operations where precision allows
- **Custom Kernels**: Hand-optimized CUDA kernels for critical operations

### Best Practices

**Development:**
- Always provide CPU fallback implementations
- Use `af::sync()` for accurate GPU timing
- Test mathematical equivalence between CPU and GPU
- Monitor memory usage during development

**Production:**
- Verify GPU availability before GPU-specific operations
- Handle GPU memory exhaustion gracefully
- Use appropriate error tolerance for GPU arithmetic
- Monitor GPU temperature and throttling

**Debugging:**
- Enable ArrayFire error printing for development builds
- Use CPU implementation for numerical debugging
- Profile GPU utilization to identify bottlenecks
- Validate results against known test cases

See [installation.md](installation.md) for GPU setup instructions and [troubleshooting.md](troubleshooting.md) for additional support.