# GPU Optimization Test Summary

## Test Configuration Details

### Dataset Information
- **Pulsar**: PSR J1733-3716
- **Parameter File**: `tests/test_data/test.par`
- **Timing File**: `tests/test_data/test.tim`
- **Number of TOAs**: 129
- **Time Span**: 4041.06 days (11.064 years)
- **Data Range**: MJD 53962.43795749 to 58012.53020181

### Model Configuration
- **Configuration File**: `test_gpu_optimization.json`
- **Test Mode**: Enabled (`"test_mode": true`)
- **Sampler**: MultiNest with 100 live points (reduced for fast testing)

### Model Elements Used
1. **Timing Model**
   - Type: Full marginalization (`"marginalise": "all"`)
   - Parameters marginalized: 5 (RAJ, DECJ, F0, F1, derived from .par file)

2. **Power Law Red Noise**
   - **days_per_coeff**: 30.0 (default)
   - **Frequencies**: 134 (calculated as floor(4041.06/30.0))
   - **Total coefficients**: 268 (134 × 2 for sine/cosine)
   - **Parameters**:
     - amplitude: log-uniform [-18, -10]
     - spectral_index: uniform [0, 7]

3. **Power Law DM Noise**
   - **days_per_coeff**: 30.0 (default)
   - **Frequencies**: 134 (calculated as floor(4041.06/30.0))
   - **Total coefficients**: 268 (134 × 2 for sine/cosine)
   - **Parameters**:
     - amplitude: log-uniform [-18, -10]
     - spectral_index: uniform [0, 7]

4. **EFAC**
   - Type: Global
   - **Parameters**:
     - global: uniform [-1, 0.7]

5. **EQUAD**
   - Type: Global
   - **Parameters**:
     - global: log-uniform [-9, -3]

### Total Model Dimensions
- **Sampling Parameters**: 6
  - Red noise amplitude (1)
  - Red noise spectral index (1)
  - DM noise amplitude (1)
  - DM noise spectral index (1)
  - EFAC global (1)
  - EQUAD global (1)
- **Marginalized Parameters**: 5 (timing model)

## Expected Test Results

### Likelihood Test
- **Test Configuration**: Fixed cube values [0.91, 0.85, 0.85, 0.9, 0.8]
- **Expected Likelihood**: 929.961754
- **Tolerance**: 1e-6
- **Status**: Updated from legacy value (956.0983174) due to coefficient calculation improvements

### GPU Optimization Details
- **Optimization Applied**: ArrayFire broadcasting for column multiplication
- **File Modified**: `likelihoods/gpu_functions.cpp:150`
- **Change**: Replaced sequential for-loop with `af::array NT = gpu_data::total_matrix_ * afNoise;`
- **Expected Performance Gain**: 10-50x for matrix multiplication operation

## Historical Context

### Legacy vs Current Implementation
- **Legacy**: Used fixed `num_freqs` parameter (typically 33)
- **Current**: Uses time-span dependent `days_per_coeff` parameter (default 30.0)
- **Impact**: More scientifically accurate frequency sampling
- **Coefficient Calculation**: `num_freqs = floor(time_span_days / days_per_coeff)`

### Test Value Change Explanation
The likelihood test value changed from 956.0983174 to 929.961754 due to:
1. **Frequency Grid Update**: Implementation of `days_per_coeff` parameter
2. **More Coefficients**: 134 frequencies vs previous fixed count
3. **Scientific Accuracy**: Better representation of noise processes over observation time span

## GPU Hardware Details
- **Device**: NVIDIA GeForce RTX 4090
- **Memory**: 24210 MB
- **Compute Capability**: 8.9
- **ArrayFire Version**: 3.9.0
- **CUDA Version**: Runtime 12.2, Driver 535.230.02

## Reproducibility Information

### Environment Setup
```bash
export TEMPO2=/home/aparthas/software/apps/tempo2
export PATH=/home/aparthas/software/apps/tempo2/bin:$PATH
```

### Build Commands
```bash
make temponest
make temponest-install
```

### Test Command
```bash
tempo2 -gr temponest -f tests/test_data/test.par tests/test_data/test.tim -Cfile test_gpu_optimization.json
```

### Expected Output Fragments
```
Time span: 4041.06 days
Red Noise: 134 frequencies (30 days/coeff, total coeffs: 268)
DM Noise: 134 frequencies (30 days/coeff, total coeffs: 268)
GPU Acceleration: ON
Likelihood test passed
```

## Version Information
- **Date**: 2025-06-26
- **TempoNest Version**: GPU-optimized with coefficient fixes
- **Key Improvements**: 
  - Time-span dependent frequency calculation
  - GPU broadcasting optimization
  - Enhanced parameter robustness

## Future Test Considerations
- For performance benchmarking, increase live_points to 500+ for realistic comparison
- Test with larger datasets to see full GPU speedup benefits
- Compare against CPU-only version for performance validation
- Test with different `days_per_coeff` values for sensitivity analysis