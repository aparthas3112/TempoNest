# TempoNest Configuration Guide

## Overview

TempoNest uses JSON configuration files to define analysis parameters, model elements, and sampler settings. This guide covers the complete configuration system with syntax, validation, and best practices.

## Table of Contents

1. [JSON Structure](#json-structure)
2. [Global Settings](#global-settings)
3. [Sampler Configuration](#sampler-configuration)
4. [Model Elements](#model-elements)
5. [Parameter Types](#parameter-types)
6. [Validation and Error Handling](#validation-and-error-handling)

---

## JSON Structure

TempoNest configuration files have three main sections:

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
        }
      ]
    }
  ]
}
```

---

## Global Settings

The `globals` section controls overall analysis behavior:

### Basic Settings

```json
{
  "globals": {
    "use_original_errors": true,
    "num_tempo2_its": 1,
    "test_mode": false
  }
}
```

**Parameters:**

- **`use_original_errors`** (boolean, default: `true`)
  - `true`: Use original TOA uncertainties from .tim file
  - `false`: Use Tempo2-fitted error estimates
  - **Recommendation**: `true` for most analyses

- **`num_tempo2_its`** (integer, default: `1`)
  - Number of Tempo2 fitting iterations before analysis
  - Range: 0-5 (higher values rarely needed)
  - **Recommendation**: `1` for standard analysis, `0` for debugging

- **`test_mode`** (boolean, default: `false`)
  - Enables additional debugging output and validation
  - Useful for development and troubleshooting
  - **Recommendation**: `false` for production runs

---

## Sampler Configuration

The `sampler` section controls the nested sampling process:

### MultiNest Settings

```json
{
  "sampler": {
    "type": "multinest",
    "sample": true,
    "live_points": 1000,
    "efficiency": 0.1,
    "evidence_tolerance": 0.5,
    "max_iterations": 0,
    "output_root": "results/TNest-",
    "resume": false,
    "verbose": true
  }
}
```

**Core Parameters:**

- **`type`** (string, required): `"multinest"`
  - Currently the only supported sampler
  - Future: HMC sampling support planned

- **`sample`** (boolean, default: `true`)
  - Whether to run sampling or just setup
  - Set to `false` for configuration validation

- **`live_points`** (integer, default: `1000`)
  - Number of live points for nested sampling
  - **Testing**: 500 (20-60 minutes)
  - **Standard**: 1000 (1-2 hours)
  - **Publication**: 4000 (4-8 hours)

- **`efficiency`** (float, default: `0.1`)
  - MultiNest efficiency parameter
  - **Fast**: 0.3 (less exploration)
  - **Balanced**: 0.1 (recommended)
  - **Thorough**: 0.01 (slower but more complete)

### Output Settings

- **`output_root`** (string, default: `"TNest-"`)
  - Prefix for all output files
  - Can include directory path: `"results/analysis1-"`
  - Generated files: `{root}post_equal_weights.dat`, `{root}stats.dat`, etc.

- **`resume`** (boolean, default: `false`)
  - Resume from previous run using `{root}resume.dat`
  - Useful for interrupted long runs

- **`verbose`** (boolean, default: `true`)
  - Control output verbosity
  - Set to `false` for batch processing

### Advanced Parameters

- **`evidence_tolerance`** (float, default: `0.5`)
  - Stopping criterion for evidence calculation
  - Smaller values = more thorough but slower

- **`max_iterations`** (integer, default: `0`)
  - Maximum sampling iterations (0 = no limit)
  - Safety stop for runaway analyses

---

## Model Elements

The `elements` array defines the physical and instrumental models:

### Element Structure

```json
{
  "name": "Element Name",
  "element_specific_options": "value",
  "parameters": [
    {
      "name": "parameter_name",
      "prior_type": "uniform",
      "min_value": 0.0,
      "max_value": 1.0,
      "include": true,
      "description": "Parameter description"
    }
  ]
}
```

### Supported Elements

See [models.md](models.md) for detailed information on each element:

- **Timing Model**: `"Timing Model"`
- **Noise Models**: `"Power Law Red Noise"`, `"Power Law DM Noise"`
- **Error Models**: `"EFAC"`, `"EQUAD"`, `"ECORR"`
- **Solar Wind**: `"Deterministic Solar Wind"`, `"Stochastic Solar Wind"`

### Element-Specific Options

**Power Law Models:**
```json
{
  "name": "Power Law Red Noise",
  "days_per_coeff": 30.0,
  "parameters": [...]
}
```

**ECORR Models:**
```json
{
  "name": "ECORR",
  "parameters": [
    {
      "name": "per_backend",
      "flag": "-B",
      "epoch_window": 0.00011574,
      "min_toas_per_epoch": 1,
      "..."
    }
  ]
}
```

**Timing Model:**
```json
{
  "name": "Timing Model",
  "marginalise": "all",  // "all", "manual", false
  "fit_all": true,       // if marginalise is false
  "sigma_multiplier": 10, // if fit_all is true
  "parameters": [...]
}
```

---

## Parameter Types

### Prior Types

**Uniform Prior:**
```json
{
  "name": "parameter_name",
  "prior_type": "uniform",
  "min_value": 0.0,
  "max_value": 7.0
}
```

**Log-Uniform Prior:**
```json
{
  "name": "amplitude",
  "prior_type": "log_uniform",
  "min_value": -18,
  "max_value": -10
}
```

### Parameter Control

**Include/Exclude Parameters:**
```json
{
  "name": "spectral_index",
  "include": true,  // false to exclude from sampling
  "prior_type": "uniform",
  "min_value": 0,
  "max_value": 7
}
```

**Flag-Based Parameters:**
```json
{
  "name": "per_flag",
  "flag": "-group",  // Apply to specific flag type
  "prior_type": "uniform",
  "min_value": -1,
  "max_value": 0.7
}
```

### Common Parameter Names

**Noise Amplitudes:**
- `amplitude` - Power law noise amplitude (log₁₀)
- `log_amplitude` - Stochastic solar wind amplitude

**Spectral Properties:**
- `spectral_index` - Power law spectral index (γ)

**Error Scaling:**
- `global` - Global EFAC parameter
- `per_flag` - Per-flag EFAC/EQUAD/ECORR parameters

**Solar Wind:**
- `electron_density` - Deterministic solar wind scaling

**Timing Parameters:**
- Use Tempo2 parameter names (e.g., `F0`, `F1`, `RAJ`, `DECJ`)

---

## Validation and Error Handling

### Configuration Validation

TempoNest performs comprehensive validation:

1. **JSON Syntax**: Valid JSON format
2. **Required Fields**: All mandatory parameters present
3. **Parameter Ranges**: Sensible min/max values
4. **Element Compatibility**: Valid element combinations
5. **Flag Validation**: Flags exist in timing data

### Common Validation Errors

**Missing Required Fields:**
```json
// Error: Missing prior_type
{
  "name": "amplitude",
  "min_value": -18,
  "max_value": -10
  // "prior_type": "log_uniform"  <-- Required!
}
```

**Invalid Parameter Ranges:**
```json
// Error: min_value > max_value
{
  "name": "spectral_index",
  "prior_type": "uniform",
  "min_value": 7,     // Should be < max_value
  "max_value": 0
}
```

**Unknown Element Names:**
```json
// Error: Typo in element name
{
  "name": "Power Law Red Nois",  // Missing 'e'
  "parameters": [...]
}
```

### Best Practices

**File Organization:**
```bash
# Use descriptive filenames
config_red_noise_analysis.json
config_solar_wind_study.json
config_comprehensive_ecorr.json

# Include version or date
config_v2.json
config_2025_analysis.json
```

**Parameter Documentation:**
```json
{
  "name": "amplitude",
  "description": "Red noise amplitude in log10(s^2/Hz)",
  "prior_type": "log_uniform",
  "min_value": -18,
  "max_value": -10
}
```

**Modular Configurations:**
- Start with simple models and add complexity gradually
- Test each element individually before combining
- Use consistent parameter ranges across analyses

### Debugging Tips

**Validation Mode:**
```json
{
  "sampler": {
    "sample": false,  // Skip sampling for validation
    "verbose": true
  },
  "globals": {
    "test_mode": true  // Enable debug output
  }
}
```

**Parameter Testing:**
```json
// Test single parameters first
{
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "all",
      "parameters": []
    }
    // Add other elements one by one
  ]
}
```

## Configuration Examples

See [examples.md](examples.md) for complete, working configuration files for common analysis scenarios.

## GUI Configuration Tool

For interactive configuration creation, use the Streamlit-based GUI:

```bash
streamlit run src/python/temponest/TempoNest_JSON.py
```

See [python_tools.md](python_tools.md) for detailed GUI usage instructions.