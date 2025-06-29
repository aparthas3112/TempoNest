# TempoNest Configuration Examples

## Overview

This document provides complete, working configuration examples for common pulsar timing analysis scenarios. Each example includes the full JSON configuration and explanations of the chosen parameters.

## Table of Contents

1. [Basic Red Noise Analysis](#basic-red-noise-analysis)
2. [Comprehensive Noise Analysis](#comprehensive-noise-analysis)
3. [Solar Wind Analysis](#solar-wind-analysis)
4. [ECORR Multi-Flag Analysis](#ecorr-multi-flag-analysis)
5. [Timing Parameter Exploration](#timing-parameter-exploration)
6. [High-Precision Publication Analysis](#high-precision-publication-analysis)

---

## Basic Red Noise Analysis

A minimal configuration for detecting red noise in pulsar timing data.

### Use Case
- First-time analysis of a pulsar dataset
- Quick assessment of timing noise properties
- Testing data quality and analysis pipeline

### Configuration

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
    "live_points": 500,
    "efficiency": 0.1,
    "output_root": "basic_red_noise/TNest-",
    "resume": false,
    "verbose": true
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
          "description": "Red noise amplitude in log10(s^2/Hz)",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Red noise spectral index (gamma)",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Error scaling factor per observing group",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-group"
        }
      ]
    }
  ]
}
```

### Key Features
- **Fast execution**: 500 live points for ~30-60 minute runtime
- **Timing model marginalization**: Automatic handling of all .par file parameters
- **Per-group EFAC**: Separate error scaling for different observing setups
- **Conservative priors**: Wide ranges suitable for discovery

---

## Comprehensive Noise Analysis

Complete noise characterization including red noise, DM noise, and white noise components.

### Use Case
- Detailed noise characterization for well-observed pulsars
- Preparation for gravitational wave analyses
- Understanding systematic effects in timing data

### Configuration

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
    "output_root": "comprehensive_analysis/TNest-",
    "resume": false,
    "verbose": true
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
          "description": "Red noise amplitude in log10(s^2/Hz)",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Red noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "Power Law DM Noise",
      "days_per_coeff": 30.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "DM noise amplitude in log10(s^2/Hz)",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "DM noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Error scaling per backend group",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-group"
        }
      ]
    },
    {
      "name": "EQUAD",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Quadrature noise per backend group",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -9,
          "max_value": -3,
          "flag": "-group"
        }
      ]
    }
  ]
}
```

### Key Features
- **Dual noise processes**: Both red and DM noise components
- **Complete error model**: EFAC and EQUAD for systematic effects
- **Standard runtime**: 1000 live points for ~1-2 hour analysis
- **Publication quality**: Suitable for detailed scientific studies

---

## Solar Wind Analysis

Specialized configuration for studying solar wind effects on pulsar timing.

### Use Case
- Investigating solar wind variations in timing data
- Correcting for interplanetary medium effects
- Solar activity correlation studies

### Configuration

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
    "output_root": "solar_wind_analysis/TNest-",
    "resume": false,
    "verbose": true
  },
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "all",
      "parameters": []
    },
    {
      "name": "Deterministic Solar Wind",
      "parameters": [
        {
          "name": "electron_density",
          "description": "Solar wind electron density scaling factor",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0.0,
          "max_value": 10.0
        }
      ]
    },
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 30.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "Residual red noise after solar wind correction",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Residual red noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Error scaling per observing group",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-group"
        }
      ]
    }
  ]
}
```

### Requirements
- Tempo2 data must include `tdis2` and `ne_sw` values
- Use pulsars with significant solar wind delays (low ecliptic latitude)

### Alternative: Stochastic Solar Wind

Replace "Deterministic Solar Wind" with:

```json
{
  "name": "Stochastic Solar Wind",
  "parameters": [
    {
      "name": "log_amplitude",
      "description": "Stochastic solar wind noise amplitude",
      "prior_type": "log_uniform",
      "include": true,
      "min_value": -18.0,
      "max_value": -10.0
    }
  ]
}
```

---

## ECORR Multi-Flag Analysis

Advanced epoch correlation analysis with multiple instrumental flag types.

### Use Case
- Detailed instrumental systematics characterization
- Multi-backend/multi-frontend timing analyses
- Legacy TNECORR parameter migration

### Configuration

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
    "output_root": "ecorr_analysis/TNest-",
    "resume": false,
    "verbose": true
  },
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "all",
      "parameters": []
    },
    {
      "name": "ECORR",
      "parameters": [
        {
          "name": "per_backend",
          "description": "Backend epoch correlations",
          "prior_type": "log_uniform",
          "include": true,
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
          "description": "Frontend epoch correlations",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -9,
          "max_value": -3,
          "flag": "-fe",
          "epoch_window": 0.00011574,
          "min_toas_per_epoch": 1
        }
      ]
    },
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 30.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "Red noise amplitude after ECORR correction",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Red noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Error scaling per group",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-group"
        }
      ]
    }
  ]
}
```

### Key Features
- **Multi-flag ECORR**: Separate correlations for backends (`-B`) and frontends (`-fe`)
- **Legacy compatibility**: 10-second epoch windows matching legacy TempoNest
- **Automatic backend detection**: System discovers backend names from data
- **Shared priors**: Same prior ranges across all flags for consistency

---

## Timing Parameter Exploration

Configuration for exploring correlations between timing parameters and noise.

### Use Case
- Investigating timing parameter uncertainties
- Parameter correlation studies
- Cases where timing model is poorly constrained

### Configuration

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
    "live_points": 1500,
    "efficiency": 0.05,
    "output_root": "timing_exploration/TNest-",
    "resume": false,
    "verbose": true
  },
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "manual",
      "parameters": [
        {
          "name": "F0",
          "description": "Spin frequency",
          "prior_type": "uniform",
          "include": true,
          "min_value": -5,
          "max_value": 5
        },
        {
          "name": "F1",
          "description": "Spin frequency derivative",
          "prior_type": "uniform",
          "include": true,
          "min_value": -3,
          "max_value": 3
        },
        {
          "name": "RAJ",
          "description": "Right ascension",
          "prior_type": "uniform",
          "include": true,
          "min_value": -5,
          "max_value": 5
        },
        {
          "name": "DECJ",
          "description": "Declination",
          "prior_type": "uniform",
          "include": true,
          "min_value": -5,
          "max_value": 5
        }
      ]
    },
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 30.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "Red noise amplitude",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Red noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-group"
        }
      ]
    }
  ]
}
```

### Key Features
- **Manual timing parameter selection**: Specific parameters included in sampling
- **Higher live points**: 1500 for better parameter correlation mapping
- **Lower efficiency**: 0.05 for more thorough exploration
- **Parameter ranges**: ±5σ around Tempo2 fitted values (in units of parameter uncertainty)

---

## High-Precision Publication Analysis

Production-quality configuration for final scientific results.

### Use Case
- Final analysis for publication
- Maximum precision and thoroughness
- Complete systematic uncertainty characterization

### Configuration

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
    "live_points": 4000,
    "efficiency": 0.01,
    "evidence_tolerance": 0.1,
    "output_root": "publication_analysis/TNest-",
    "resume": true,
    "verbose": true
  },
  "elements": [
    {
      "name": "Timing Model",
      "marginalise": "all",
      "parameters": []
    },
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 20.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "Red noise amplitude in log10(s^2/Hz)",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "Red noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "Power Law DM Noise",
      "days_per_coeff": 20.0,
      "parameters": [
        {
          "name": "amplitude",
          "description": "DM noise amplitude in log10(s^2/Hz)",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -18,
          "max_value": -10
        },
        {
          "name": "spectral_index",
          "description": "DM noise spectral index",
          "prior_type": "uniform",
          "include": true,
          "min_value": 0,
          "max_value": 7
        }
      ]
    },
    {
      "name": "ECORR",
      "parameters": [
        {
          "name": "per_backend",
          "description": "Backend epoch correlations",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -9,
          "max_value": -3,
          "flag": "-B",
          "epoch_window": 0.00011574,
          "min_toas_per_epoch": 1
        }
      ]
    },
    {
      "name": "EFAC",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Error scaling per backend",
          "prior_type": "uniform",
          "include": true,
          "min_value": -1,
          "max_value": 0.7,
          "flag": "-B"
        }
      ]
    },
    {
      "name": "EQUAD",
      "parameters": [
        {
          "name": "per_flag",
          "description": "Quadrature noise per backend",
          "prior_type": "log_uniform",
          "include": true,
          "min_value": -9,
          "max_value": -3,
          "flag": "-B"
        }
      ]
    }
  ]
}
```

### Key Features
- **Maximum precision**: 4000 live points for ~4-8 hour runtime
- **Thorough exploration**: 0.01 efficiency for complete posterior mapping
- **High frequency resolution**: 20 days per coefficient for finer spectral detail
- **Complete noise model**: All major systematic effects included
- **Resume capability**: Essential for long runs that may be interrupted
- **Tight evidence tolerance**: 0.1 for maximum precision

### Expected Runtime
- **CPU**: 6-12 hours on 100-200 cores
- **GPU**: 4-8 hours on single high-end GPU (RTX 4090 class)

---

## Usage Guidelines

### Running Configurations

```bash
# Basic usage
tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json

# With MPI (recommended for large analyses)
mpirun -np 8 tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json
```

### Configuration Testing

Before long runs, test configurations with reduced settings:

1. **Validation run**: Set `"sample": false` to check configuration
2. **Quick test**: Use 100-200 live points for rapid verification
3. **Parameter check**: Use `"test_mode": true` for debugging

### Output Files

Each configuration produces:
- `{output_root}post_equal_weights.dat` - Posterior samples
- `{output_root}stats.dat` - Evidence and parameter statistics  
- `{output_root}summary.txt` - Human-readable summary
- `{output_root}live.points` - Live points during sampling

See [python_tools.md](python_tools.md) for analysis and plotting tools.