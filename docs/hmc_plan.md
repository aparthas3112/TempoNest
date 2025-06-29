# HMC Sampling Implementation Plan for TempoNest

## Overview
Implement GPU-accelerated Hamiltonian Monte Carlo (HMC) sampling using JAX/NumPyro while maintaining full compatibility with existing MultiNest workflow and achieving 10-100x speedup.

**Target Deployment**: Portable installation across GPU machines and HPC clusters with minimal dependencies.

## Phase 1: Foundation & Setup (2-3 days)

### 1.1 Dependency Installation & Testing
**Primary Goal**: Establish portable, cluster-ready installation procedure

#### CUDA Compatibility Strategy
- **Current Approach**: Proceed with existing ArrayFire CUDA version unless conflicts arise
- **Upgrade Path**: If needed, reinstall ArrayFire with recent CUDA (11.8/12.x) for JAX compatibility
- **HPC Consideration**: Use containerized deployment (Docker/Singularity) for cluster portability

#### Installation Procedure (HPC-Ready)
```bash
# Option 1: Conda-based (recommended for HPC)
conda create -n temponest-hmc python=3.10
conda activate temponest-hmc
conda install -c conda-forge jax jaxlib[cuda] numpyro arviz

# Option 2: pip with CUDA compatibility check
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install numpyro arviz optax pybind11

# Option 3: Containerized deployment (for HPC clusters)
# Dockerfile with all dependencies pre-built
```

#### Compatibility Testing Suite
- Test JAX-ArrayFire coexistence with current CUDA setup
- Validate GPU memory sharing between frameworks
- Create automated compatibility check script for new installations
- Document fallback procedures for different HPC environments

### 1.2 Sampler Architecture Abstraction
**Goal**: Clean separation between sampling algorithms with unified interface

#### Core Architecture Changes
```cpp
// New sampler factory system
class sampler_factory_t {
public:
    enum class type_t { MULTINEST, HMC_NUTS };
    static std::unique_ptr<sampler_t> create(const json_node_t& config);
};

// Abstract sampler interface (existing)
class sampler_t {
public:
    virtual void run(std::shared_ptr<model_t> model) = 0;
    virtual void output_results() = 0;
    virtual bool validate_configuration() = 0;
};
```

#### JSON Configuration Extension
```json
{
  "sampler": {
    "type": "hmc",  // or "multinest" 
    "hmc_specific": {
      "num_chains": 4,
      "num_samples": 2000,
      "num_warmup": 1000,
      "step_size": 0.01,        // Initial step size (auto-tuned)
      "target_accept_prob": 0.8, // Target acceptance probability
      "max_tree_depth": 10,     // NUTS tree depth limit
      "num_leapfrog_steps": 50  // Alternative to NUTS
    },
    "multinest_specific": {
      "live_points": 500,
      "efficiency": 0.1
      // ... existing parameters
    }
  }
}
```

#### Backward Compatibility
- Existing JSON configs without "sampler.type" default to MultiNest
- All current MultiNest parameters remain functional
- Gradual migration path for existing users

### 1.3 Evidence Computation Research & Implementation
**Goal**: Achieve evidence accuracy within 0.5 log units of MultiNest

#### Bridge Sampling Implementation
```python
# Primary evidence computation method
from numpyro.infer import MCMC, NUTS
import arviz as az
from scipy.special import logsumexp

class evidence_computer_t:
    def bridge_sampling_evidence(self, mcmc_samples, model_fn):
        """
        Compute log evidence using bridge sampling
        Target accuracy: ±0.5 log units vs MultiNest
        """
        pass
    
    def cross_validate_evidence(self, multinest_evidence, hmc_evidence):
        """Validation against MultiNest results"""
        pass
```

#### Evidence Validation Strategy
- Compare HMC bridge sampling vs MultiNest evidence on known test cases
- Implement multiple evidence estimators (WAIC, LOO, TI) for cross-validation
- Create evidence benchmark suite for different model complexities

## Phase 2: HMC Sampler Integration (3-4 days)

### 2.1 HMC Sampler Shell
**Goal**: Complete HMC sampler with identical interface to MultiNest

#### Core Implementation
```cpp
class hmc_sampler_t : public sampler_t {
private:
    // Python integration
    py::object numpyro_model_;
    py::object mcmc_runner_;
    py::object evidence_computer_;
    
    // Configuration
    hmc_settings_t settings_;
    
public:
    // Same interface as multinest_sampler_t
    void run(std::shared_ptr<model_t> model) override;
    void output_results() override;
    bool validate_configuration() override;
    
    // HMC-specific
    void initialize_jax_environment();
    void validate_gpu_availability();
    py::object create_numpyro_model(const model_space_t& model_space);
};
```

#### Robust Error Handling
```cpp
class hmc_validator_t {
public:
    void check_gpu_availability() {
        // Fail fast if GPU not available for HMC
        if (!jax_gpu_available()) {
            throw std::runtime_error(
                "HMC sampler requires GPU but none available. "
                "Use 'multinest' sampler for CPU-only analysis."
            );
        }
    }
    
    void validate_jax_installation() {
        // Check JAX/NumPyro installation
        // Provide helpful error messages for missing dependencies
    }
};
```

### 2.2 JAX Likelihood Bridge
**Goal**: Seamless C++ ↔ JAX integration with automatic differentiation

#### Critical Implementation Details
```cpp
// Python-C++ bridge for likelihood evaluation
class jax_likelihood_bridge_t {
private:
    std::shared_ptr<temponest_v1_t> cpp_likelihood_;
    
public:
    // Create JAX-compatible likelihood function
    py::object create_jax_likelihood() {
        return py::cpp_function([this](py::array_t<double> params) {
            // Convert numpy array to std::vector
            std::vector<double> cpp_params = numpy_to_vector(params);
            
            // Call C++ likelihood
            double result = (*cpp_likelihood_)(model_space_, cpp_params);
            
            return result;
        });
    }
    
    // Validation: ensure C++ and JAX give identical results
    void validate_consistency() {
        // Test multiple parameter sets
        // Require machine precision agreement (< 1e-12)
    }
};
```

#### Handling JAX Differentiability Issues
```python
# Replace problematic conditional operations
def make_jax_compatible(original_likelihood):
    """Convert conditional logic to JAX-friendly operations"""
    
    # Example: Replace if/else with jnp.where
    # Original: if totCoeff > 0: update_diagonal()
    # JAX-friendly:
    diagonal_update = jnp.where(
        totCoeff > 0,
        powercoeff_inverse,
        jnp.zeros_like(powercoeff_inverse)
    )
    
    return jax_likelihood
```

### 2.3 Model Translation System
**Goal**: Perfect preservation of TempoNest model definitions in NumPyro

#### Model Element Translation
```python
class model_translator_t:
    """Convert TempoNest JSON models to NumPyro format"""
    
    def translate_timing_model(self, timing_config):
        """Convert timing model parameters to NumPyro priors"""
        pass
    
    def translate_noise_models(self, noise_configs):
        """Convert red noise, DM noise to NumPyro distributions"""
        # Preserve exact power law specifications
        # red_amp: log-uniform(-18, -10) → dist.Uniform(-18, -10)
        # red_gamma: uniform(0, 7) → dist.Uniform(0, 7)
        pass
    
    def translate_efac_equad(self, error_configs):
        """Convert EFAC/EQUAD models with per-flag support"""
        pass
    
    def translate_solar_wind(self, sw_configs):
        """Convert deterministic/stochastic solar wind models"""
        pass
```

#### Parameter Bounds & Transformations
```python
# Ensure identical parameter handling
def create_numpyro_model(temponest_config):
    """Generate NumPyro model matching TempoNest exactly"""
    
    # Example: Red noise amplitude (log-uniform)
    red_amp_log = numpyro.sample(
        "red_amplitude_log", 
        dist.Uniform(low=-18.0, high=-10.0)
    )
    red_amp = jnp.exp(red_amp_log)  # Transform to linear space
    
    # Likelihood evaluation
    likelihood_value = jax_likelihood_bridge(red_amp, ...)
    numpyro.factor("likelihood", likelihood_value)
```

## Phase 3: Results & Validation (2-3 days)

### 3.1 Results Conversion Pipeline
**Goal**: Seamless integration with existing TempoNest output ecosystem

#### Chain Processing
```cpp
class results_converter_t {
public:
    // Convert NumPyro MCMC samples to TempoNest format
    std::vector<parameter_stats_t> process_chains(py::object mcmc_samples) {
        // Extract parameter chains
        // Compute means, standard deviations, correlations
        // Generate MAP and maximum likelihood estimates
        // Calculate effective sample sizes and R-hat diagnostics
    }
    
    // Generate standard TempoNest outputs
    void write_output_files() {
        // chains.txt (compatible with existing corner plot tools)
        // summary.txt (parameter statistics)
        // evidence.txt (bridge sampling results)
        // diagnostics.txt (convergence information)
    }
    
    // Integration with existing visualization
    void generate_corner_plots() {
        // Use existing plot_corner.py
        // Ensure chain format compatibility
    }
};
```

#### Output Format Preservation
```bash
# Maintain existing file structure
chains/
├── chains.txt           # MCMC samples (TempoNest format)
├── summary.txt          # Parameter statistics  
├── evidence.txt         # Log evidence estimate
├── diagnostics.txt      # HMC-specific diagnostics
└── corner_plot.pdf      # Generated by existing tools
```

### 3.2 Comprehensive Validation Suite
**Goal**: Rigorous testing ensuring scientific accuracy

#### Mathematical Consistency Tests
```cpp
class validation_suite_t {
public:
    // Test 1: Likelihood evaluation consistency
    void test_likelihood_consistency() {
        // Generate test parameter sets
        // Compare C++ vs JAX likelihood (require < 1e-12 difference)
        // Test edge cases (boundary parameters, complex models)
    }
    
    // Test 2: Prior sampling validation  
    void test_prior_distributions() {
        // Sample from NumPyro priors vs TempoNest priors
        // Statistical tests for distribution agreement
    }
    
    // Test 3: Parameter transformation consistency
    void test_parameter_bounds() {
        // Verify log-uniform, uniform transformations identical
        // Test boundary behavior
    }
};
```

#### Statistical Validation Protocol
```python
def cross_validate_posteriors(multinest_results, hmc_results):
    """
    Comprehensive posterior comparison
    Target: Parameter means within 1%, correlations within 0.05
    """
    
    # Kolmogorov-Smirnov tests for each parameter
    for param in parameters:
        ks_stat, p_value = stats.ks_2samp(
            multinest_results[param], 
            hmc_results[param]
        )
        assert p_value > 0.05, f"Parameter {param} distributions differ"
    
    # Parameter mean comparison (within 1 sigma)
    for param in parameters:
        mn_mean, mn_std = multinest_results[param].mean(), multinest_results[param].std()
        hmc_mean, hmc_std = hmc_results[param].mean(), hmc_results[param].std()
        
        diff = abs(mn_mean - hmc_mean)
        combined_std = np.sqrt(mn_std**2 + hmc_std**2)
        assert diff < combined_std, f"Parameter {param} means differ significantly"
```

#### Evidence Validation
```python
def validate_evidence_computation():
    """Test evidence accuracy vs MultiNest"""
    
    test_cases = [
        "simple_red_noise",      # 2 parameters
        "red_dm_efac",          # 5 parameters  
        "complex_solar_wind",   # 10+ parameters
    ]
    
    for case in test_cases:
        multinest_evidence = run_multinest_evidence(case)
        hmc_evidence = run_hmc_bridge_sampling(case)
        
        evidence_diff = abs(multinest_evidence - hmc_evidence)
        assert evidence_diff < 0.5, f"Evidence differs by {evidence_diff} log units"
```

### 3.3 Continuous Testing Integration
**Goal**: Automated validation preventing regressions

#### Test Suite Architecture
```bash
# Automated testing after significant changes
./run_hmc_validation.sh

# Test components:
1. Installation compatibility check
2. GPU availability validation  
3. Mathematical consistency tests
4. Statistical posterior validation
5. Evidence computation accuracy
6. Performance benchmarking
7. Output format compatibility
```

#### Performance Regression Testing
```cpp
class performance_monitor_t {
public:
    void benchmark_speedup() {
        // Measure HMC vs MultiNest runtime
        // Target: 10-100x speedup for realistic problems
        // Track performance regressions
    }
    
    void monitor_memory_usage() {
        // GPU memory utilization
        // CPU-GPU transfer overhead
        // Memory leak detection
    }
};
```

## Technical Implementation Details

### Installation & Deployment Strategy

#### HPC Cluster Compatibility
```bash
# Containerized deployment for maximum portability
FROM nvidia/cuda:12.1-devel-ubuntu22.04

# Install TempoNest dependencies
RUN apt-get update && apt-get install -y \
    build-essential cmake gfortran \
    libgsl-dev libfftw3-dev

# Install Python environment
RUN conda create -n temponest python=3.10
RUN conda install -c conda-forge jax jaxlib[cuda] numpyro

# Build TempoNest with HMC support
COPY . /temponest
RUN cd /temponest && ./configure --enable-hmc && make
```

#### Module System Integration
```bash
# Environment module for HPC systems
module load cuda/12.1
module load python/3.10
module load temponest/hmc-gpu

# Or Spack package manager integration
spack install temponest+hmc+gpu
```

### Memory Management Strategy

#### GPU Memory Coordination
```cpp
class gpu_memory_manager_t {
private:
    size_t arrayfire_pool_size_;
    size_t jax_pool_size_;
    
public:
    void initialize_memory_pools() {
        // Allocate separate pools for ArrayFire and JAX
        // Avoid memory fragmentation
        // Monitor usage and adjust dynamically
    }
    
    void synchronize_contexts() {
        // Coordinate between ArrayFire likelihood and JAX sampling
        // Minimize CPU-GPU transfers
    }
};
```

#### Data Transfer Optimization
```cpp
// Minimize Python-C++ data transfers
class efficient_data_bridge_t {
public:
    // Batch parameter evaluations
    py::array evaluate_likelihood_batch(py::array param_batch) {
        // Process multiple parameter sets at once
        // Reduce transfer overhead
    }
    
    // In-place operations where possible
    void update_parameters_inplace(py::array params) {
        // Avoid unnecessary copying
    }
};
```

## Risk Mitigation & Contingency Plans

### Technical Risk Management

#### CUDA Version Conflicts
```bash
# Contingency plan for CUDA incompatibility
if [ cuda_version_conflict ]; then
    echo "Reinstalling ArrayFire with compatible CUDA version..."
    # Automated ArrayFire reinstallation procedure
    reinstall_arrayfire_with_cuda12()
fi
```

#### JAX Installation Issues
```python
# Robust installation checking
def validate_jax_installation():
    try:
        import jax
        import jax.numpy as jnp
        import numpyro
        
        # Test GPU availability
        assert jax.devices('gpu'), "No GPU devices found"
        
        # Test basic operations
        test_array = jnp.array([1.0, 2.0, 3.0])
        assert jnp.sum(test_array) == 6.0
        
        print("JAX installation validated successfully")
        
    except Exception as e:
        raise RuntimeError(f"JAX installation issue: {e}")
```

### Scientific Risk Management

#### Evidence Computation Fallbacks
```python
# Multiple evidence estimation methods
class robust_evidence_computer_t:
    def compute_evidence(self, samples):
        try:
            # Primary: Bridge sampling
            return self.bridge_sampling_evidence(samples)
        except Exception:
            # Fallback: WAIC approximation
            return self.waic_evidence(samples)
```

#### Convergence Monitoring
```python
def monitor_hmc_convergence():
    """Robust convergence diagnostics"""
    
    # R-hat statistic (should be < 1.01)
    rhat = az.rhat(mcmc_samples)
    
    # Effective sample size (should be > 400)
    ess = az.ess(mcmc_samples)
    
    # Divergent transitions (should be < 1%)
    divergences = mcmc_samples.get_diagnostics()['divergences']
    
    if any([rhat > 1.01, ess < 400, divergences > 0.01]):
        raise Warning("HMC convergence issues detected")
```

## Success Criteria & Deliverables

### Quantitative Success Metrics
1. **Mathematical Consistency**: C++ and JAX likelihood evaluations agree to machine precision (< 1e-12)
2. **Statistical Accuracy**: HMC posteriors match MultiNest within 1% for parameter means, 0.05 for correlations
3. **Evidence Reliability**: HMC evidence estimates within 0.5 log units of MultiNest
4. **Performance Target**: 10-100x speedup over MultiNest on realistic problems (10+ parameters)
5. **Scientific Validation**: Accurate recovery of injected signals in simulated datasets
6. **Convergence Reliability**: R-hat < 1.01, ESS > 400 for all parameters

### Deliverables
- ✅ Complete HMC sampler implementation with GPU acceleration
- ✅ Bridge sampling evidence computation with validation
- ✅ Comprehensive automated testing suite
- ✅ HPC-ready installation and deployment procedures
- ✅ Performance benchmarking and comparison studies
- ✅ Migration guide and documentation for existing users
- ✅ Containerized deployment for cluster environments

### User Experience Goals
- **Seamless Integration**: Existing JSON configs work with minimal changes
- **Clear Error Messages**: Helpful guidance when GPU/dependencies unavailable
- **Identical Output**: Same file formats and visualization tools work unchanged
- **Easy Installation**: One-command setup for common HPC environments
- **Robust Fallbacks**: Graceful degradation when components unavailable

## Timeline & Milestones

**Week 1: Foundation**
- Day 1-2: Dependency installation and compatibility testing
- Day 3-4: Sampler architecture abstraction
- Day 5-7: Evidence computation implementation

**Week 2: Integration** 
- Day 1-3: JAX likelihood bridge and model translation
- Day 4-5: Results conversion and output compatibility
- Day 6-7: Comprehensive validation and testing

**Ongoing**: Performance optimization and documentation

This timeline is flexible based on technical challenges and validation requirements. Priority is ensuring scientific accuracy and robustness over speed of implementation.