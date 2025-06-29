# TempoNest Documentation

## Overview
TempoNest is a Bayesian pulsar timing analysis tool with GPU acceleration, implementing time-span dependent noise modeling and advanced parameter handling for robust scientific analysis.

## Table of Contents
1. [Installation and Setup](#installation-and-setup)
2. [Configuration System](#configuration-system)
3. [Model Elements](#model-elements)
4. [Technical Architecture](#technical-architecture)
5. [Performance and Optimization](#performance-and-optimization)
6. [Tools and Utilities](#tools-and-utilities)
7. [Examples and Use Cases](#examples-and-use-cases)
8. [Development Guidelines](#development-guidelines)
9. [Recent Updates](#recent-updates)
10. [Troubleshooting](#troubleshooting)

## Installation and Setup

### Dependencies
Before building, ensure these are installed and environment variables are set:
- `$TEMPO2` - Path to Tempo2 installation (required)
- MultiNest library with Fortran bindings
- GSL, BLAS/LAPACK, Eigen3, MPI compiler (mpicc/mpicxx)
- Optional: ArrayFire for GPU acceleration

### Build Commands
```bash
# Initial setup (run once after cloning)
./autogen.sh
./configure

# Standard build
make temponest

# Install plugin to Tempo2
make temponest-install

# Run tests
make test

# Clean build
make clean
```

### Configure Options
Key configure flags:
- `--enable-debug` - Debug build with symbols
- `--with-eigen=/path` - Specify Eigen headers location
- `--with-arrayfire=/path` - Enable GPU support
- `--with-tempo2-plug-dir=/path` - Custom plugin install directory

### Testing
```bash
# Run all tests (builds first if needed)
make test

# Manual test run with Tempo2
$TEMPO2/bin/tempo2 -gr temponest -f tests/test_data/test.par tests/test_data/test.tim -cfile tests/test_data/test.json
```

## Configuration System

### JSON Structure
TempoNest uses JSON configuration files with three main sections:

```json
{
  "globals": {
    "use_original_errors": true,
    "num_tempo2_its": 1,
    "test_mode": false
  },
  "sampler": {
    "type": "multinest",
    "sample": true,
    "live_points": 1000,
    "efficiency": 0.1,
    "output_root": "results/TNest-"
  },
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "all",
      "parameters": []
    },
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 30.0,
      "parameters": [
        {
          "name": "amplitude",
          "prior_type": "log_uniform",
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "prior_type": "uniform",
          "min_value": 0,
          "max_value": 7
        }
      ]
    }
  ]
}
```

### Global Settings
- `use_original_errors`: Use original TOA uncertainties vs Tempo2-fitted errors
- `num_tempo2_its`: Number of Tempo2 iterations before analysis
- `test_mode`: Enable additional debugging output

### Sampler Configuration
- `type`: "multinest" (primary sampler)
- `live_points`: 500 (testing), 1000 (standard), 4000 (publication)
- `efficiency`: 0.3 (fast), 0.1 (balanced), 0.01 (thorough)
- `output_root`: Output file prefix

## Model Elements

### Timing Model
Handles Tempo2 pulsar parameters with three modes:

#### 1. Full Marginalization (Recommended)
```json
{
  "name": "Timing Model",
  "marginalise": "all",
  "parameters": []
}
```
Automatically marginalizes over all fitted parameters from .par file.

#### 2. Fit All Parameters
```json
{
  "name": "Timing Model",
  "marginalise": false,
  "fit_all": true,
  "sigma_multiplier": 10,
  "parameters": []
}
```
Includes all timing parameters as sampled parameters with ±10σ ranges.

#### 3. Manual Parameter Selection
```json
{
  "name": "Timing Model",
  "marginalise": "manual",
  "parameters": [
    {
      "name": "F0",
      "prior_type": "uniform",
      "include": true,
      "min_value": -5,
      "max_value": 5
    }
  ]
}
```

### Power Law Red Noise
Models timing correlated noise with power-law spectrum:

**Mathematical Model**: S(f) = A² × (f/f₁ᵧᵣ)^(-γ) / (12π²)

```json
{
  "name": "Power Law Red Noise",
  "days_per_coeff": 30.0,
  "parameters": [
    {
      "name": "amplitude",
      "description": "log₁₀ amplitude of power law noise",
      "prior_type": "log_uniform",
      "min_value": -18,
      "max_value": -10
    },
    {
      "name": "spectral_index",
      "description": "spectral index (γ, red noise)",
      "prior_type": "uniform",
      "min_value": 0,
      "max_value": 7
    }
  ]
}
```

**Key Parameters**:
- `days_per_coeff`: Time span per frequency coefficient (default: 30.0 days)
- `amplitude`: Log₁₀ amplitude in s²/Hz
- `spectral_index`: Power law index (γ > 0 for red noise)

### Power Law DM Noise
Models dispersion measure variations:

**Mathematical Model**: S(f) = A² × f₁ᵧᵣ^(-3) × (f × 365.25)^(-γ) / T_span

```json
{
  "name": "Power Law DM Noise",
  "days_per_coeff": 30.0,
  "parameters": [
    {
      "name": "amplitude",
      "description": "log₁₀ amplitude of DM noise",
      "prior_type": "log_uniform",
      "min_value": -18,
      "max_value": -10
    },
    {
      "name": "spectral_index",
      "description": "spectral index (γ, DM noise)",
      "prior_type": "uniform",
      "min_value": 0,
      "max_value": 7
    }
  ]
}
```

### EFAC (Error Factor)
Scales TOA uncertainties by constant factor:

#### Global EFAC
```json
{
  "name": "EFAC",
  "parameters": [
    {
      "name": "global",
      "prior_type": "uniform",
      "min_value": -1,
      "max_value": 0.7
    }
  ]
}
```

#### Per-Flag EFAC
```json
{
  "name": "EFAC",
  "parameters": [
    {
      "name": "per_flag",
      "prior_type": "uniform",
      "min_value": -1,
      "max_value": 0.7,
      "flag": "-group"
    }
  ]
}
```

**Mathematical Model**: σ'ᵢ = σᵢ × 10^(EFAC)

### EQUAD (Quadrature Noise)
Adds white noise in quadrature with TOA uncertainties:

```json
{
  "name": "EQUAD",
  "parameters": [
    {
      "name": "per_flag",
      "prior_type": "log_uniform",
      "min_value": -9,
      "max_value": -3,
      "flag": "-group"
    }
  ]
}
```

**Mathematical Model**: σ'ᵢ = √(σᵢ² + 10^(2×EQUAD))

### ECORR (Epoch Correlations)
Models timing noise correlations within observation epochs:

#### Multi-Flag Configuration
```json
{
  "name": "ECORR",
  "parameters": [
    {
      "name": "per_backend",
      "prior_type": "log_uniform",
      "min_value": -9,
      "max_value": -3,
      "flag": "-B",
      "epoch_window": 0.00011574,
      "min_toas_per_epoch": 1
    }
  ]
}
```

#### Advanced ECORR Features
- **Multi-flag Support**: Separate ECORR elements for different flag types (`-B`, `-fe`, `-f`, `-sys`, `-chan`, `-group`)
- **Automatic Backend Detection**: System discovers backends within each flag from data
- **Legacy Compatibility**: 10-second epoch windows, minimum 1 TOA per epoch
- **Shared Priors**: Same prior ranges applied across all selected flags

**Mathematical Model**: 
- **Variance**: σ² = 10^(2×log₁₀(ECORR)) where ECORR is in microseconds
- **Covariance**: C = σ² U U^T where U is the quantization matrix

#### Legacy TNECORR Migration
**From Legacy .par entries**:
```
TNECORR -B 10CM 0.0257802
TNECORR -B uwl_10CM 0.000232653
TNECORR -fe UWL 0.222656
TNECORR -f KAT_MKBF 0.0622733
```

**To Current JSON**: Select multiple flags in GUI with shared priors - system auto-detects backends and creates equivalent parameter structure.

### Solar Wind Models

#### Deterministic Solar Wind
Corrects timing residuals for systematic solar wind variations:

```json
{
  "name": "Deterministic Solar Wind",
  "parameters": [
    {
      "name": "electron_density",
      "description": "Solar wind electron density scaling",
      "prior_type": "uniform",
      "min_value": 0.0,
      "max_value": 10.0
    }
  ]
}
```

**Mathematical Model**: δt = (ne - ne_ref) × tdis2

#### Stochastic Solar Wind
Adds frequency-dependent white noise scaled by solar wind delay:

```json
{
  "name": "Stochastic Solar Wind",
  "parameters": [
    {
      "name": "log_amplitude",
      "description": "Log₁₀ amplitude of stochastic solar wind noise",
      "prior_type": "log_uniform",
      "min_value": -18.0,
      "max_value": -10.0
    }
  ]
}
```

**Usage Notes**:
- Typically use deterministic OR stochastic, not both
- Requires Tempo2 to provide `tdis2` and `ne_sw` values
- Full GUI support with validation and helpful tooltips

## Technical Architecture

### Core Design Pattern
TempoNest uses a **compositional model system** where analysis models are built by combining independent model elements:

```
model_t (wrapper)
├── model_space_t (container for elements)
│   ├── timing_model_t (Tempo2 parameters)
│   ├── efac_t (error scaling)
│   ├── equad_t (quadrature noise)
│   ├── ecorr_t (epoch correlations)
│   └── power_law_*_noise_t (stochastic processes)
└── likelihood_t (computation engine)
    └── temponest_v1_t (Gaussian likelihood)
```

### Key Components

#### Model Elements (`model_elements/`)
Each element handles a specific aspect of pulsar timing:
- Inherit from `model_element_t` base class
- Manage their own parameters via `parameter_t` objects
- Implement `update_residuals()` to modify Tempo2 data
- Use JSON configuration for setup

#### Parameter System (`types/parameter.h`)
- Supports uniform, log-uniform priors
- Automatic indexing for sampler parameter vectors
- Parent-child relationships with model elements
- Prior bounds and transformations

#### Samplers (`samplers/`)
- Factory pattern creates samplers from JSON config
- `multinest_t` implements nested sampling
- Settings classes handle configuration and validation
- Abstract interface allows adding new sampling methods

#### Plugin Integration
TempoNest integrates with Tempo2 as a dynamically loaded plugin:
- Entry point: `graphicalInterface()` in `TempoNestPlugin.cpp`
- Uses Tempo2's pulsar data structures directly
- Installs as `.t2` shared library in Tempo2 plugins directory
- Command line: `tempo2 -gr temponest -f par.file tim.file -Cfile config.json`

### GPU Implementation

#### ArrayFire Integration
All GPU operations use ArrayFire for optimal performance:
```cpp
// Example: Broadcasting operation
af::array result = matrix * vector;  // Automatic broadcasting
af::array chol_L;
af::cholesky(chol_L, TNT);           // GPU Cholesky decomposition
```

#### Memory Management
Static GPU data structure prevents repeated allocation:
```cpp
namespace gpu_data {
    static af::array total_matrix_;
    static bool initialized_ = false;
}
```

#### GPU Function Development
```cpp
#ifdef HAVE_ARRAYFIRE
// GPU implementation
af::array gpu_result = optimized_gpu_function(gpu_data);
#else
// CPU fallback
Eigen::VectorXd cpu_result = standard_cpu_function(cpu_data);
#endif
```

## Performance and Optimization

### GPU Acceleration
- **Hardware Tested**: NVIDIA GeForce RTX 4090 (24GB VRAM)
- **ArrayFire Version**: 3.9.0 with CUDA backend
- **Typical Speedup**: 1-2x overall system performance vs 300+ CPU cores
- **Memory Bandwidth**: >1000 GB/s vs ~100 GB/s per CPU socket

### Scaling Behavior
- **Parameter Dimensions**: Efficient for >20 parameter problems
- **Time Complexity**: O(N³) for Cholesky operations, where N = observations + coefficients
- **Memory Usage**: Proportional to N² for covariance matrices

### MultiNest Performance
**Runtime scaling**: O(nlive × log(evidence_range))
- nlive=500: 20-60 min (testing)
- nlive=1000: 1-2 hours (standard)  
- nlive=4000: 4-8 hours (publication)

**Efficiency settings**:
- 0.3: Fast convergence, less exploration
- 0.1: Balanced (recommended)
- 0.01: Thorough, slower

**Typical acceptance rates**: 5-15%

### Recent GPU Optimizations (2025-06)

#### Major GPU Performance Improvements
**Problem**: Sequential column multiplication in likelihood calculation was destroying GPU parallelization.

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
- Expected 10-50x speedup for matrix operations
- GPU utilization improved from ~1% to near-100% for this operation
- Verified mathematically identical results to CPU implementation

#### Coefficient Calculation Overhaul
**Problem**: Legacy implementation used fixed frequency counts (`num_freqs = 33`) independent of observation time span, creating scientifically inconsistent noise models.

**Solution**: Implemented time-span dependent `days_per_coeff` parameter:
```cpp
// New implementation
int num_freqs = static_cast<int>(std::floor(time_span_days / days_per_coeff));
```

**Impact**:
- Default `days_per_coeff = 30.0` provides scientifically accurate frequency sampling
- Example: 4041-day span → 134 frequencies (vs fixed 33 previously)
- Total coefficients: 268 red + 268 DM = 536 vs ~132 previously
- Mathematical consistency with legacy CPU versions restored

## Tools and Utilities

### Configuration GUI (`TempoNest_JSON.py`)
Streamlit-based interface for JSON configuration generation:

#### Key Features
- **Three timing model modes**: Marginalise (default), Fit All, Manual
- **Enhanced solar wind support** with proper default ranges
- **ECORR multi-flag support**: Select multiple flags with shared priors
- **Improved defaults**: Per-flag EFAC/EQUAD, 500 live points, local directory output
- **Element type validation**: Parameter ranges reset correctly when changing element types

#### Usage
```bash
streamlit run TempoNest_JSON.py
```

### Enhanced Corner Plot Tool (`plot_corner.py`)

#### Basic Usage
```bash
python plot_corner.py chains.txt  # Auto-parameter detection with LaTeX labels
```

#### Advanced Comparison Mode
Overlay posteriors from multiple analysis results:

```bash
# Compare two analyses with same parameter structure
python plot_corner.py chains1.txt --compare /path/to/dir2/ \
    --param-map-1 red_amp,red_gamma,efac \
    --param-map-2 red_amp,red_gamma,efac

# Compare with different parameter orderings
python plot_corner.py chains1.txt --compare /path/to/dir2/ \
    --param-map-1 0,1,2 \
    --param-map-2 2,0,1

# Multiple comparisons with custom labels
python plot_corner.py main_chains.txt --compare dir2/ --compare dir3/ \
    --param-map-1 red_amp,red_gamma,efac \
    --param-map-2 red_amp,red_gamma,efac --param-map-2 0,1,2 \
    --label-1 "GPU Analysis" --label-2 "CPU Test" --label-2 "Legacy"
```

#### Key Features
- **Multiple Dataset Support**: Compare up to 10 datasets simultaneously
- **Flexible Parameter Mapping**: Handle different parameter orders and naming conventions
- **Enhanced File Detection**: Automatically finds parameter name files
- **Custom Labels**: User-defined legend labels
- **Visual Distinction**: 10-color palette with transparency
- **Parameter Recognition**: Enhanced with solar wind parameter recognition and LaTeX formatting

## Examples and Use Cases

### Basic Red Noise Analysis
```json
{
  "globals": {"use_original_errors": true, "num_tempo2_its": 1},
  "sampler": {"type": "multinest", "sample": true, "live_points": 500, "efficiency": 0.1},
  "elements": [
    {"name": "Timing Model", "marginalise": "all", "parameters": []},
    {"name": "Power Law Red Noise", "parameters": [
      {"name": "amplitude", "prior_type": "log_uniform", "include": true, "min_value": -18, "max_value": -10},
      {"name": "spectral_index", "prior_type": "uniform", "include": true, "min_value": 0, "max_value": 7}
    ]},
    {"name": "EFAC", "parameters": [
      {"name": "per_flag", "prior_type": "uniform", "include": true, "min_value": -1, "max_value": 0.7, "flag": "-group"}
    ]}
  ]
}
```

### Solar Wind Analysis
```json
{
  "globals": {"use_original_errors": true, "num_tempo2_its": 1},
  "sampler": {"type": "multinest", "sample": true, "live_points": 500, "efficiency": 0.1},
  "elements": [
    {"name": "Timing Model", "marginalise": "all", "parameters": []},
    {"name": "Deterministic Solar Wind", "parameters": [
      {"name": "electron_density", "prior_type": "uniform", "include": true, "min_value": 0.0, "max_value": 10.0}
    ]},
    {"name": "EFAC", "parameters": [
      {"name": "per_flag", "prior_type": "uniform", "include": true, "min_value": -1, "max_value": 0.7, "flag": "-group"}
    ]}
  ]
}
```

### Comprehensive ECORR Analysis
```json
{
  "globals": {"use_original_errors": true, "num_tempo2_its": 1},
  "sampler": {"type": "multinest", "sample": true, "live_points": 1000, "efficiency": 0.1},
  "elements": [
    {"name": "Timing Model", "marginalise": "all", "parameters": []},
    {"name": "ECORR", "parameters": [
      {"name": "per_backend", "prior_type": "log_uniform", "include": true, "min_value": -9, "max_value": -3, 
       "flag": "-B", "epoch_window": 0.00011574, "min_toas_per_epoch": 1}
    ]},
    {"name": "ECORR", "parameters": [
      {"name": "per_backend", "prior_type": "log_uniform", "include": true, "min_value": -9, "max_value": -3,
       "flag": "-fe", "epoch_window": 0.00011574, "min_toas_per_epoch": 1}
    ]},
    {"name": "Power Law Red Noise", "parameters": [
      {"name": "amplitude", "prior_type": "log_uniform", "include": true, "min_value": -18, "max_value": -10},
      {"name": "spectral_index", "prior_type": "uniform", "include": true, "min_value": 0, "max_value": 7}
    ]},
    {"name": "EFAC", "parameters": [
      {"name": "per_flag", "prior_type": "uniform", "include": true, "min_value": -1, "max_value": 0.7, "flag": "-group"}
    ]}
  ]
}
```

## Development Guidelines

### Adding Model Elements
1. Inherit from `model_element_t`
2. Implement `update_residuals()` method
3. Register in factory pattern (`model_space.cpp`)
4. Add GPU implementation if applicable (`gpu_functions.cpp`)

### Testing Protocol
1. **Unit Tests**: Verify mathematical equivalence between CPU/GPU
2. **Integration Tests**: Full likelihood calculations with known datasets
3. **Performance Tests**: Benchmark against baseline implementations
4. **Regression Tests**: Ensure coefficient calculations remain consistent

### File Organization

#### Core Implementation
- `likelihoods/gpu_functions.cpp`: GPU-accelerated likelihood calculations
- `model_elements/`: All noise and timing model implementations
- `tests/test_likelihood.cpp`: Regression testing framework

#### Configuration and Tools
- `TempoNest_JSON.py`: Streamlit GUI for configuration generation
- `plot_corner.py`: Enhanced corner plot tool with comparison capabilities
- `GPU_OPTIMIZATION_TEST_SUMMARY.md`: Test documentation and expected results

## Recent Updates (2024-2025)

### ECORR Implementation and Build System Fixes

#### Critical Bug Fix: Undefined Symbol Error
**Problem Resolved**: Fixed `undefined symbol: _ZTI7ecorr_t` error when loading TempoNest plugin
- **Root Cause**: `ecorr.cpp` was missing from `Makefile.am` build sources 
- **Solution**: Added `model_elements/ecorr.cpp` to `common_sources` in build system
- **Impact**: ECORR functionality now fully available

#### ECORR Model Implementation (Complete)
**Status**: Full implementation of epoch correlations with legacy compatibility

**Core Features**:
- **Mathematical Model**: σ² = 10^(2×log₁₀(ECORR)) where ECORR is in microseconds
- **Epoch Detection**: Configurable time windows (default: 10 seconds, legacy compatible)
- **Backend Discovery**: Automatic detection of backends within each flag type
- **Quantization Matrix**: Binary 0/1 matrices for epoch assignment (U[i,j] = 1 if TOA i belongs to epoch j)

**Legacy Compatibility**:
- **Epoch Window**: 10 seconds (matches legacy TempoNest)
- **Minimum TOAs**: 1 TOA per epoch (matches legacy behavior)
- **Mathematical Equivalence**: Identical covariance model C = σ² U U^T

### JUMP Parameter Analysis (Comprehensive Review)
**Status**: Complete verification of scientific equivalence between legacy and current implementations

**Mathematical Consistency**:
- **Legacy doJumpMargin**: User-controlled via `doJumpMargin=0/1` parameter
- **Current Implementation**: Always analytically marginalized (automatic)
- **Scientific Results**: Mathematically identical (SVD-based marginalization preserved exactly)

**Robustness Improvements**:
- **Parameter Scanner**: Modern C++ container-based approach replaces `fitinfo` dependency
- **Error Prevention**: Comprehensive bounds checking prevents crashes on complex .par files
- **JUMP Types Supported**: JUMP, FDJUMP (frequency-dependent), FDJUMPDM (DM-like)

### Mathematical Model Verification (2024)
**Status**: Complete verification of mathematical consistency between current and legacy implementations

**Models Verified**:
- **Red Noise**: Identical power-law spectrum `S(f) = A² × (f/f₁ᵧᵣ)^(-γ) / (12π²)`
- **DM Noise**: Identical `S(f) = A² × f₁ᵧᵣ^(-3) × (f × 365.25)^(-γ) / T_span`  
- **Solar Wind**: Identical deterministic `δt = (ne - ne_ref) × tdis2` and stochastic formulations
- **Scaling Factors**: All constants (`f1yr = 1/3.16e7`, normalization factors) preserved exactly

**Conclusion**: No mathematical differences between versions - scientific results fully consistent

## Troubleshooting

### Common Issues

#### Build Problems
- **Missing TEMPO2**: Ensure `$TEMPO2` environment variable is set
- **ArrayFire not found**: Use `--with-arrayfire=/path` or disable GPU with `--without-arrayfire`
- **Eigen3 headers**: Usually auto-detected, specify with `--with-eigen=/path` if needed

#### Runtime Errors
- **Plugin not found**: Run `make temponest-install` to install to Tempo2 plugins directory
- **Segmentation faults**: Often due to malformed .par files - use debug build with `--enable-debug`
- **JSON parsing errors**: Validate JSON syntax, check required fields

#### Performance Issues
- **Slow convergence**: Increase live points or adjust efficiency parameter
- **Memory usage**: Reduce problem size or increase system RAM
- **GPU not utilized**: Verify ArrayFire installation and CUDA drivers

### Reproducibility and Testing

#### Reference Test Case
- **Dataset**: PSR J1733-3716 (129 TOAs, 4041-day span)
- **Configuration**: `test_gpu_optimization.json`
- **Expected Likelihood**: 929.961754 ± 1e-6
- **Frequency Count**: 134 red + 134 DM = 268 total frequencies

#### Environment Requirements
```bash
export TEMPO2=/path/to/tempo2
export PATH=$TEMPO2/bin:$PATH
make temponest && make temponest-install
tempo2 -gr temponest -f test.par test.tim -Cfile config.json
```

### Known Limitations

#### Current Limitations
1. **Diagonal Updates**: Sequential GPU implementation pending optimization
2. **Memory Pooling**: Static allocation may be suboptimal for varying problem sizes
3. **Batch Processing**: Single dataset processing limits GPU utilization

#### Future Optimization Opportunities
1. **Phase 2**: Diagonal vectorization (expected 3-5x additional speedup)
2. **Phase 3**: Memory pooling and batch processing (expected 2-3x additional speedup)
3. **Advanced**: Multi-GPU support for extremely large problems

---

This documentation reflects the state as of June 2025 with major coefficient, GPU optimizations, and ECORR implementation completed and verified.