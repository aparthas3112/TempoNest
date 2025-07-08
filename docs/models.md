# TempoNest Model Elements

## Overview

TempoNest uses a compositional model system where pulsar timing analyses are built by combining independent model elements. Each element handles a specific aspect of pulsar timing physics, from timing model parameters to various noise processes.

## Table of Contents

1. [Timing Model](#timing-model)
2. [Noise Models](#noise-models)
   - [Power Law Red Noise](#power-law-red-noise)
   - [Power Law DM Noise](#power-law-dm-noise)
3. [Error Models](#error-models)
   - [EFAC (Error Factor)](#efac-error-factor)
   - [EQUAD (Quadrature Noise)](#equad-quadrature-noise)
   - [ECORR (Epoch Correlations)](#ecorr-epoch-correlations)
4. [Solar Wind Models](#solar-wind-models)
   - [Deterministic Solar Wind](#deterministic-solar-wind)
   - [Stochastic Solar Wind](#stochastic-solar-wind)

---

## Timing Model

The timing model handles Tempo2 pulsar parameters with three operational modes for maximum flexibility.

### 1. Full Marginalization (Recommended)

Automatically marginalizes over all fitted parameters from the .par file:

```json
{
  "name": "Timing Model",
  "marginalise": "all",
  "parameters": []
}
```

**Use Case**: Standard analysis where timing parameters are well-constrained by data.

### 2. Fit All Parameters

Includes all timing parameters as sampled parameters with automatic range determination:

```json
{
  "name": "Timing Model",
  "marginalise": false,
  "fit_all": true,
  "sigma_multiplier": 10,
  "parameters": []
}
```

**Parameters**:
- `sigma_multiplier`: Range factor (±Nσ) around Tempo2 fitted values (default: 10)

**Use Case**: When timing parameters are poorly constrained or when exploring parameter correlations.

### 3. Manual Parameter Selection

Fine-grained control over which parameters to include:

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
    },
    {
      "name": "F1",
      "prior_type": "uniform",
      "include": true,
      "min_value": -3,
      "max_value": 3
    }
  ]
}
```

**Use Case**: Custom analysis requiring specific parameter subset or non-standard prior ranges.

---

## Noise Models

### Power Law Red Noise

Models timing-correlated noise with a power-law spectrum, commonly observed in pulsar timing data.

**Mathematical Model**: 
```
S(f) = A² × (f/f₁ᵧᵣ)^(-γ) / (12π²)
```

Where:
- `A`: Amplitude in s²/Hz
- `γ`: Spectral index (γ > 0 for red noise)
- `f₁ᵧᵣ = 1/(1 year) = 1/3.16×10⁷ s⁻¹`: Reference frequency

**Configuration**:
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
  - Determines frequency sampling: `num_freqs = floor(time_span_days / days_per_coeff)`
  - Example: 4041-day span → 134 frequencies
- `amplitude`: Log₁₀ amplitude, typically in range [-18, -10]
- `spectral_index`: Power law index, typical red noise has γ ∈ [1, 7]

### Power Law DM Noise

Models dispersion measure (DM) variations with power-law spectrum.

**Mathematical Model**:
```
S(f) = A² × f₁ᵧᵣ^(-3) × (f × 365.25)^(-γ) / T_span
```

Where the `f₁ᵧᵣ^(-3)` term accounts for DM scaling and `365.25` converts to annual units.

**Configuration**:
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

**Physical Interpretation**: DM variations can arise from changes in the interstellar medium electron density along the line of sight.

---

## Error Models

### EFAC (Error Factor)

Scales Time-of-Arrival (TOA) uncertainties by a constant multiplicative factor to account for systematic underestimation of measurement errors.

**Mathematical Model**: 
```
σ'ᵢ = σᵢ × 10^(EFAC)
```

#### Global EFAC

Applies single scaling factor to all TOAs:

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

Separate scaling factors for different instrumental setups:

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

**Common Flags**: `-group`, `-sys`, `-fe` (frontend), `-be` (backend)

### EQUAD (Quadrature Noise)

Adds white noise in quadrature with TOA uncertainties to account for unmodeled systematic effects.

**Mathematical Model**: 
```
σ'ᵢ = √(σᵢ² + 10^(2×EQUAD))
```

**Configuration**:
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

**Units**: EQUAD values are log₁₀(microseconds)

### ECORR (Epoch Correlations)

Models timing noise correlations within observation epochs, addressing instrumental systematics that affect multiple TOAs collected during the same observing session.

**Mathematical Model**: 
- **Variance**: `σ² = 10^(2×log₁₀(ECORR))` where ECORR is in microseconds
- **Covariance Matrix**: `C = σ² U U^T` where U is the quantization matrix
- **Quantization Matrix**: Binary matrix where `U[i,j] = 1` if TOA i belongs to epoch j

#### Basic ECORR Configuration

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

#### Advanced Multi-Flag ECORR

Handle multiple instrumental configurations:

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
},
{
  "name": "ECORR",
  "parameters": [
    {
      "name": "per_backend",
      "prior_type": "log_uniform",
      "min_value": -9,
      "max_value": -3,
      "flag": "-fe",
      "epoch_window": 0.00011574,
      "min_toas_per_epoch": 1
    }
  ]
}
```

**Key Features**:
- **Multi-flag Support**: Separate ECORR elements for different flag types (`-B`, `-fe`, `-f`, `-sys`, `-chan`, `-group`)
- **Automatic Backend Detection**: System discovers backends within each flag from data
- **Legacy Compatibility**: 10-second epoch windows (0.00011574 days), minimum 1 TOA per epoch
- **Shared Priors**: Same prior ranges applied across all selected flags

**Supported Flags**:
- `-B`: Backend systems
- `-fe`: Frontend systems  
- `-f`: Frequency/filter specifications
- `-sys`: System identifiers
- `-chan`: Channel information
- `-group`: Arbitrary groupings

#### Legacy TNECORR Migration

**From Legacy .par entries**:
```
TNECORR -B 10CM 0.0257802
TNECORR -B uwl_10CM 0.000232653
TNECORR -fe UWL 0.222656
TNECORR -f KAT_MKBF 0.0622733
```

**To Current JSON**: Use multi-flag ECORR configuration - system auto-detects backends and creates equivalent parameter structure.

---

## Solar Wind Models

### Deterministic Solar Wind

Corrects timing residuals for systematic solar wind variations using Tempo2-provided electron density measurements.

**Mathematical Model**: 
```
δt = (ne - ne_ref) × tdis2
```

Where:
- `ne`: Solar wind electron density from Tempo2
- `ne_ref`: Reference electron density  
- `tdis2`: Solar wind delay coefficient from Tempo2

**Configuration**:
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

**Requirements**: Tempo2 must provide `tdis2` and `ne_sw` values in timing data.

### Stochastic Solar Wind

Models unmodeled solar wind variations using either white noise or Gaussian process models. Both approaches use the ν⁻² frequency scaling characteristic of dispersion measure variations.

**Two Operating Modes:**

#### White Noise Mode

Adds frequency-dependent white noise scaled by solar wind delay.

**Mathematical Model**: 
```
Noise variance ∝ A × (tdis2/ne_sw)²
```

**Configuration**:
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

#### Gaussian Process Mode

Models solar wind variations as correlated noise with power-law temporal spectrum.

**Mathematical Model**:
```
P(f_temporal) = A² × f_temporal^(-γ)
Frequency scaling: ν⁻² (fixed, from tdis2/ne_sw)
```

**Configuration**:
```json
{
  "name": "Stochastic Solar Wind",
  "days_per_coeff": 30.0,
  "parameters": [
    {
      "name": "log10_A_sw",
      "description": "Log₁₀ amplitude of stochastic solar wind GP noise",
      "prior_type": "log_uniform",
      "min_value": -18.0,
      "max_value": -10.0
    },
    {
      "name": "gamma_sw",
      "description": "Spectral index for stochastic solar wind GP",
      "prior_type": "uniform",
      "min_value": 0.0,
      "max_value": 7.0
    }
  ]
}
```

**Key Parameters for GP Mode**:
- `days_per_coeff`: Time span per frequency coefficient (default: 30.0 days)
- `log10_A_sw`: Log₁₀ amplitude, typically in range [-18, -10]
- `gamma_sw`: Spectral index, typical values γ ∈ [0, 7]

**Physical Interpretation**: 
- **White Mode**: Accounts for uncorrelated solar wind fluctuations
- **GP Mode**: Models temporally correlated solar wind variations (e.g., solar cycle effects)
- **Combined Mode**: Both modes can be used simultaneously for comprehensive modeling

**Enterprise Compatibility**: 
- Frequency scaling (ν⁻²) matches Enterprise's solar wind implementation
- Uses Tempo2's pre-computed solar wind geometry (tdis2 parameter)
- Temporal power spectrum follows standard pulsar timing conventions

**Usage Notes**:
- Can be used alongside deterministic solar wind for comprehensive modeling
- Requires Tempo2 to provide `tdis2` and `ne_sw` values
- GP mode provides more sophisticated modeling of correlated solar wind variations
- Amplitude units are log₁₀(s²) for both modes

---

## Model Combination Guidelines

### Typical Analysis Configurations

**Basic Noise Analysis**:
- Timing Model (marginalized)
- Power Law Red Noise
- EFAC (per-flag)
- EQUAD (per-flag)

**Comprehensive Analysis**:
- Timing Model (marginalized)  
- Power Law Red Noise
- Power Law DM Noise
- ECORR (multi-flag)
- EFAC (per-flag)
- EQUAD (per-flag)

**Solar Wind Study**:
- Timing Model (marginalized)
- Deterministic Solar Wind
- Power Law Red Noise
- EFAC (per-flag)

### Parameter Scaling Guidelines

**Time Span Dependencies**:
- `days_per_coeff`: Adjust based on data span and desired frequency resolution
- Longer spans → more coefficients → higher computational cost
- Typical values: 20-50 days depending on data characteristics

**Prior Ranges**:
- **Amplitude priors**: Adjust based on pulsar and noise expectations
- **Spectral indices**: Red noise typically γ ∈ [1, 7], DM noise similar
- **Error scaling**: EFAC typically near 1 (log₁₀ ≈ 0), EQUAD/ECORR depend on instrumentation

See [examples.md](examples.md) for complete configuration examples and [configuration.md](configuration.md) for JSON syntax details.