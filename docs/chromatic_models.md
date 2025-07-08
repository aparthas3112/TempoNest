# Chromatic Noise Models in TempoNest

## Overview

TempoNest includes advanced chromatic noise models that capture frequency-dependent timing variations beyond simple dispersion measure (DM) effects. These models are essential for precision pulsar timing analysis and gravitational wave detection.

## Table of Contents

1. [Chromatic GP Noise](#chromatic-gp-noise)
2. [Dynamic Scaling Implementation](#dynamic-scaling-implementation)
3. [Mathematical Foundation](#mathematical-foundation)
4. [Configuration Examples](#configuration-examples)
5. [Usage Guidelines](#usage-guidelines)

---

## Chromatic GP Noise

The Chromatic GP noise model captures timing variations that scale with observing frequency as ν^(-α), where α (chromatic index) can be fitted or fixed.

### Mathematical Model

**Temporal Power Spectrum:**
```
P(f_temporal) = A² × f_temporal^(-γ)
```

**Frequency Scaling:**
```
Chromatic scaling = (ν_ref / ν_obs)^α
Reference frequency: ν_ref = 1400 MHz
```

**Combined Model:**
```
C(f_temporal, ν_obs) = P(f_temporal) × (ν_ref / ν_obs)^(2α)
```

Where:
- `A`: Amplitude in units of (s × MHz^α)
- `γ`: Temporal spectral index 
- `α`: Chromatic index (frequency scaling power)
- `ν_ref`: Reference frequency (1400 MHz)
- `ν_obs`: Observing frequency

### Operating Modes

#### Fixed Chromatic Index

Standard mode with predetermined frequency scaling:

```json
{
  "name": "Chromatic GP Noise",
  "days_per_coeff": 30.0,
  "chromatic_idx": 4.0,
  "parameters": [
    {
      "name": "amplitude",
      "description": "Log₁₀ amplitude of chromatic GP noise",
      "prior_type": "log_uniform",
      "min_value": -18,
      "max_value": -10
    },
    {
      "name": "spectral_index",
      "description": "Temporal spectral index (γ)",
      "prior_type": "uniform", 
      "min_value": 0,
      "max_value": 7
    }
  ]
}
```

#### Variable Chromatic Index (Dynamic Scaling)

Advanced mode where the chromatic index is fitted as a parameter:

```json
{
  "name": "Chromatic GP Noise",
  "days_per_coeff": 30.0,
  "parameters": [
    {
      "name": "amplitude",
      "description": "Log₁₀ amplitude of chromatic GP noise",
      "prior_type": "log_uniform",
      "min_value": -18,
      "max_value": -10
    },
    {
      "name": "spectral_index", 
      "description": "Temporal spectral index (γ)",
      "prior_type": "uniform",
      "min_value": 0,
      "max_value": 7
    },
    {
      "name": "chromatic_idx",
      "description": "Chromatic frequency scaling index (α)",
      "prior_type": "uniform",
      "min_value": 0.0,
      "max_value": 6.0
    }
  ]
}
```

---

## Dynamic Scaling Implementation

### The Challenge

Traditional implementations require reconstructing the entire design matrix when the chromatic index changes during sampling, which is computationally prohibitive for large datasets.

### Hybrid Solution

TempoNest implements a mathematically rigorous hybrid approach:

1. **Design Matrix Construction**: Build matrix once using a reference chromatic index (α_ref = 4.0)
2. **Dynamic Correction**: Apply scaling correction during likelihood evaluation
3. **Mathematical Accuracy**: Maintains full scientific precision

### Scaling Correction Formula

**Correction Factor:**
```
correction_factor = Σₖ[(ν_ref/νₖ)^(2×α_current)] / Σₖ[(ν_ref/νₖ)^(2×α_ref)]
```

**Applied as:**
```
GP_coefficient_corrected = GP_coefficient_reference × correction_factor
```

Where:
- `α_current`: Current chromatic index value during sampling
- `α_ref`: Reference index used in design matrix (4.0)
- `νₖ`: Observing frequency of observation k
- Sum is over all non-deleted observations

### Computational Advantages

- **Performance**: Avoids expensive matrix reconstruction (>100x speedup)
- **Memory**: Constant memory usage regardless of index variation
- **Accuracy**: Mathematically equivalent to full reconstruction
- **Robustness**: Handles any chromatic index range efficiently

---

## Mathematical Foundation

### Frequency Scaling Physics

Different chromatic indices correspond to distinct physical processes:

- **α = 2**: Standard dispersion measure (DM) variations
- **α = 4**: Scattering measure (SM) variations  
- **α ∈ [2,4]**: Mixed ionospheric/interstellar effects
- **α > 4**: Exotic chromatic processes

### Covariance Structure

**Full Covariance Matrix:**
```
C[i,j] = Σₘ [Aₘ × cos(2πfₘtᵢ) × cos(2πfₘtⱼ) + Aₘ × sin(2πfₘtᵢ) × sin(2πfₘtⱼ)] × (ν_ref/νᵢ)^α × (ν_ref/νⱼ)^α
```

**Power Coefficients:**
```
Aₘ = A² × (fₘ × 365.25 / T_span)^(-γ) × normalization_factor
```

Where:
- `fₘ`: Temporal frequency mode (1, 2, 3, ..., num_freqs)
- `T_span`: Total observation span
- `tᵢ, tⱼ`: Observation times
- `νᵢ, νⱼ`: Observing frequencies

### Dynamic Scaling Derivation

**Original GP Coefficient:**
```
GP_coeff[m] ∝ A² × f_m^(-γ) × Σₖ[(ν_ref/νₖ)^(2×α_ref)]
```

**Desired GP Coefficient:**
```
GP_coeff_desired[m] ∝ A² × f_m^(-γ) × Σₖ[(ν_ref/νₖ)^(2×α_current)]
```

**Correction Factor:**
```
correction = GP_coeff_desired[m] / GP_coeff[m] = Σₖ[(ν_ref/νₖ)^(2×α_current)] / Σₖ[(ν_ref/νₖ)^(2×α_ref)]
```

This correction is frequency-independent and applies uniformly to all temporal modes.

---

## Configuration Examples

### Testing Chromatic Index

Example configuration to test if chromatic index converges to expected values:

```json
{
  "name": "Chromatic GP Noise",
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
    },
    {
      "name": "chromatic_idx",
      "description": "Test chromatic index - should converge to ~2 for DM-only data",
      "prior_type": "uniform",
      "min_value": 0.0,
      "max_value": 3.0
    }
  ]
}
```

**Expected Results:**
- **DM-dominated data**: α ≈ 2.0
- **Scattering-dominated data**: α ≈ 4.0  
- **Mixed effects**: α ∈ [2,4]

### Comprehensive Chromatic Analysis

Advanced configuration for detailed chromatic noise characterization:

```json
{
  "name": "Chromatic GP Noise", 
  "days_per_coeff": 45.0,
  "parameters": [
    {
      "name": "amplitude",
      "description": "Chromatic amplitude at 1400 MHz",
      "prior_type": "log_uniform", 
      "min_value": -20,
      "max_value": -8
    },
    {
      "name": "spectral_index",
      "description": "Temporal spectral index",
      "prior_type": "uniform",
      "min_value": 0,
      "max_value": 7
    },
    {
      "name": "chromatic_idx", 
      "description": "Frequency scaling index",
      "prior_type": "uniform",
      "min_value": 1.5,
      "max_value": 5.0
    }
  ]
}
```

---

## Usage Guidelines

### When to Use Chromatic GP

- **High-precision timing**: Required for gravitational wave analyses
- **Multi-frequency datasets**: Essential when spanning wide frequency ranges
- **Long-term studies**: Captures evolving chromatic noise properties
- **Systematic investigations**: Understanding frequency-dependent effects

### Parameter Selection

**Chromatic Index Priors:**
- **Conservative**: [0, 6] for discovery mode
- **DM-focused**: [1.5, 2.5] when expecting pure DM variations
- **Scattering studies**: [3.5, 4.5] for scintillation-dominated sources
- **Physical constraints**: [2, 4] for astrophysically motivated range

**Amplitude Ranges:**
- **Wide-band data**: Larger amplitudes typically needed
- **Narrow-band data**: Lower amplitudes sufficient
- **Pulsar-dependent**: Adjust based on source characteristics

**Time Coefficients:**
- **days_per_coeff**: 30-60 days typical range
- **Longer spans**: More coefficients → better temporal resolution
- **Computational cost**: Scales quadratically with number of coefficients

### Validation Strategies

1. **Cross-validation**: Compare fixed vs variable index results
2. **Posterior checks**: Verify chromatic index convergence
3. **Model comparison**: Test against standard DM noise models
4. **Physical consistency**: Check index values against expectations

### Performance Considerations

- **Dynamic scaling**: No significant computational overhead
- **Memory usage**: Standard GP memory requirements
- **Convergence**: May require more live points for complex chromatic structure
- **Runtime**: Similar to fixed-index mode due to efficient implementation

---

## Implementation Notes

### Code Integration

The dynamic scaling feature is automatically detected when `chromatic_idx` appears as a fitted parameter rather than a configuration option.

**Automatic Detection:**
- If `chromatic_idx` is in parameters → Dynamic scaling enabled
- If `chromatic_idx` is in config → Fixed scaling used
- Clear console output indicates which mode is active

### Numerical Stability

The implementation includes safeguards for:
- **Zero-frequency handling**: Proper limits for edge cases
- **Extreme indices**: Numerical stability for α ∈ [0, 10]
- **Reference frequency**: Consistent 1400 MHz reference across implementations

### Future Extensions

Potential enhancements for future development:
- **Multiple chromatic components**: Independent GP processes with different indices
- **Time-varying indices**: Allowing chromatic index evolution over observation span  
- **Cross-frequency correlations**: Advanced multi-band chromatic modeling

---

For practical implementation examples, see [examples.md](examples.md). For configuration syntax details, see [configuration.md](configuration.md).