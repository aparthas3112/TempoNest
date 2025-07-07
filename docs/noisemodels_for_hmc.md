# Noise Model Gradients for HMC Implementation

## Overview

This document provides the mathematical framework for computing gradients of the TempoNest likelihood with respect to noise parameters for Hamiltonian Monte Carlo (HMC) sampling. The timing model parameters remain analytically marginalized as in the current implementation.

## Core Likelihood Structure

The TempoNest likelihood has the form:
```
log L = -½(tdet + jointdet + freq_det + timelike - freqlike) + log_prior_densities
```

Where:
- `tdet = -∑ log(N_i)` (white noise determinant)
- `jointdet = log|M^T N M + Λ^(-1)|` (joint determinant)  
- `freq_det = ∑ log(λ_i)` (noise coefficient determinant)
- `timelike = r^T N r` (time domain quadratic form)
- `freqlike = (M^T N r)^T (M^T N M + Λ^(-1))^(-1) (M^T N r)` (frequency domain quadratic form)

**Key matrices:**
- `N = diag(1/σ_i²)`: White noise precision matrix
- `M`: Design matrix `[U_timing | M_red | M_dm | M_ecorr]`
- `Λ = diag(λ_i)`: Noise coefficient covariance matrix
- `r`: Timing residuals

## Gradient Computation Strategy

For HMC, we need `∇_θ log L` where `θ` represents noise parameters. The gradient splits into contributions from each term:

```
∇_θ log L = -½(∇_θ tdet + ∇_θ jointdet + ∇_θ freq_det + ∇_θ timelike - ∇_θ freqlike) + ∇_θ log_priors
```

## 1. White Noise Parameters (EFAC/EQUAD)

### EFAC (Error Scaling)

**Parameter transformation:**
```cpp
σ_i² → σ_i² * EFAC²  // For observations with matching flag
N_i = 1/σ_i² → N_i / EFAC²
```

**Gradients:**

**1. Time domain contribution:**
```
∂tdet/∂EFAC = ∂(-∑ log N_i)/∂EFAC = ∂(-∑ log(N_i⁰/EFAC²))/∂EFAC = 2n_flag/EFAC
```
where `n_flag` is the number of observations with this flag.

```
∂timelike/∂EFAC = ∂(r^T N r)/∂EFAC = -2/EFAC³ ∑_{i∈flag} r_i² N_i⁰
```

**2. Joint determinant contribution:**
```
∂jointdet/∂EFAC = Tr[(M^T N M + Λ^(-1))^(-1) ∂(M^T N M)/∂EFAC]
```
where:
```
∂(M^T N M)/∂EFAC = -2/EFAC³ ∑_{i∈flag} M_i M_i^T N_i⁰
```

**3. Frequency domain contribution:**
```
∂freqlike/∂EFAC = 2(M^T N r)^T (M^T N M + Λ^(-1))^(-1) ∂(M^T N r)/∂EFAC
                  - (M^T N r)^T (M^T N M + Λ^(-1))^(-1) ∂(M^T N M)/∂EFAC (M^T N M + Λ^(-1))^(-1) (M^T N r)
```

**Prior contribution:**
```
∂log_prior/∂EFAC = -1/EFAC  // For log-uniform prior
∂log_prior/∂EFAC = 0        // For uniform prior  
```

### EQUAD (Additional White Noise)

**Parameter transformation:**
```
σ_i² → σ_i² + EQUAD²  // For observations with matching flag
N_i = 1/(σ_i⁰² + EQUAD²)
```

**Gradients:**
```
∂N_i/∂EQUAD = -2*EQUAD/(σ_i⁰² + EQUAD²)² = -2*EQUAD*N_i²

∂tdet/∂EQUAD = -∑_{i∈flag} 1/N_i * ∂N_i/∂EQUAD = 2*EQUAD*n_flag

∂timelike/∂EQUAD = ∑_{i∈flag} r_i² * ∂N_i/∂EQUAD = -2*EQUAD ∑_{i∈flag} r_i² N_i²
```

## 2. Red Noise Parameters

### Red Noise Amplitude (Log-uniform prior)

**Parameter enters via power spectrum coefficients:**
```cpp
λ_red,i = (A_red² / 12π²) * (f_i/f_1yr)^(-γ) / (T * normalization)
```

**Gradients:**
```
∂λ_red,i/∂A_red = 2*A_red/12π² * (f_i/f_1yr)^(-γ) / (T * normalization) = 2λ_red,i/A_red

∂freq_det/∂A_red = ∑_i ∂log(λ_red,i)/∂A_red = 2*n_red_freqs/A_red

∂jointdet/∂A_red = Tr[(M^T N M + Λ^(-1))^(-1) ∂Λ^(-1)/∂A_red]
```
where:
```
∂Λ^(-1)/∂A_red = -diag(∂λ_red,i/∂A_red / λ_red,i²) = -2/(A_red * λ_red,i) in red noise block
```

**Frequency domain:**
```
∂freqlike/∂A_red = -(M^T N r)^T (M^T N M + Λ^(-1))^(-1) ∂Λ^(-1)/∂A_red (M^T N M + Λ^(-1))^(-1) (M^T N r)
```

**Prior contribution:**
```
∂log_prior/∂A_red = -1/A_red  // Log-uniform prior
```

### Red Noise Spectral Index (Uniform prior)

**Parameter enters via frequency dependence:**
```cpp
λ_red,i ∝ f_i^(-γ)
```

**Gradients:**
```
∂λ_red,i/∂γ = -λ_red,i * log(f_i/f_1yr)

∂freq_det/∂γ = -∑_i log(f_i/f_1yr)

∂log_prior/∂γ = 0  // Uniform prior
```

## 3. DM Noise Parameters

DM noise follows identical mathematics to red noise, but with frequency scaling modified by dispersion measure:

```cpp
M_dm,ij = cos/sin(2π f_i t_j) / (ν_j²)  // Additional 1/ν² factor
```

Gradients are identical to red noise case.

## 4. ECORR Parameters

### ECORR Amplitudes (Uniform priors)

**Parameter enters via block-diagonal additions to Λ:**
```cpp
λ_ecorr,epoch = ECORR_amp²  // For each unique epoch
```

**Gradients:**
```
∂λ_ecorr/∂ECORR_amp = 2*ECORR_amp

∂freq_det/∂ECORR_amp = 2*n_epochs/ECORR_amp

∂log_prior/∂ECORR_amp = 0  // Uniform prior
```

## 5. Solar Wind Parameters

### Deterministic Solar Wind
Modifies residuals directly: `r → r - M_sw * α_sw`

```
∂timelike/∂α_sw = -2*(r - M_sw*α_sw)^T N M_sw
∂freqlike/∂α_sw = 2*(M^T N r)^T (M^T N M + Λ^(-1))^(-1) M_sw^T N
```

### Stochastic Solar Wind
Adds to white noise variance like EQUAD.

## Implementation Strategy

### Efficient Computation

**1. Reuse existing matrix factorizations:**
```cpp
// Already computed in likelihood:
Eigen::LLT<Eigen::MatrixXd> llt(M^T N M + Λ^(-1));
Eigen::VectorXd alpha = llt.solve(M^T N r);

// For gradients, reuse:
Eigen::MatrixXd inv_TNT = llt.solve(Eigen::MatrixXd::Identity(size, size));
```

**2. Sparse gradient updates:**
```cpp
// EFAC only affects observations with specific flags
// Red noise only affects specific frequency components
// Use sparse operations where possible
```

**3. Numerical stability:**
```cpp
// Use log-space arithmetic for very small/large values
// Monitor condition numbers of matrix systems
// Implement gradient checking for validation
```

### Validation Framework

**1. Finite difference checking:**
```cpp
double finite_diff_grad = (L(θ + ε) - L(θ - ε)) / (2*ε);
double analytical_grad = compute_gradient(θ);
assert(abs(finite_diff_grad - analytical_grad) < tolerance);
```

**2. Gradient magnitude monitoring:**
```cpp
// Flag unusually large gradients that might indicate numerical issues
// Compare gradient norms across different parameter types
```

**3. Consistency tests:**
```cpp
// Verify ∇_θ L · direction ≈ directional derivative
// Check that likelihood and gradient evaluate consistently
```

## Parameter Transformation Considerations

### From MultiNest Unit Hypercube to HMC Physical Space

**Current (MultiNest):**
```cpp
// u ∈ [0,1] → physical parameter
double red_amp = pow(10, u * (log_max - log_min) + log_min);  // Log-uniform
double red_index = u * (max_index - min_index) + min_index;   // Uniform
```

**HMC approach:**
```cpp
// Work directly in physical space: red_amp ∈ [10^-18, 10^-10]
// Add proper prior densities instead of Jacobian terms
```

### Prior Density Gradients

**Log-uniform prior:** `p(x) ∝ 1/x` for `x ∈ [x_min, x_max]`
```cpp
double log_density = -log(x) - log(log(x_max) - log(x_min));
double grad_log_density = -1.0/x;
```

**Uniform prior:** `p(x) ∝ 1` for `x ∈ [x_min, x_max]`  
```cpp
double log_density = -log(x_max - x_min);
double grad_log_density = 0.0;
```

## Expected Performance

### Gradient Computation Cost
- **EFAC/EQUAD**: O(n_obs) per parameter
- **Red/DM noise**: O(n_freq) per parameter  
- **ECORR**: O(n_epochs) per parameter
- **Total**: O(n_obs + n_freq + n_epochs) ≪ O(n_obs²) for likelihood

### Memory Requirements
- **Additional storage**: ~3× likelihood memory (gradients + intermediate matrices)
- **Matrix reuse**: Leverage existing Cholesky factorization
- **Sparse operations**: Exploit parameter locality

### Convergence Expectations
- **Effective sample size**: ~10× improvement over MultiNest for correlated parameters
- **Exploration efficiency**: Gradient guidance particularly effective for red noise correlations
- **Warmup period**: ~1000 samples for step size adaptation

This mathematical framework provides the foundation for implementing robust, efficient HMC sampling of TempoNest noise parameters while preserving the sophisticated timing model marginalization that makes TempoNest scientifically powerful.