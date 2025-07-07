# CLAUDE.md

TempoNest is a Bayesian pulsar timing analysis tool that integrates with Tempo2 as a plugin. Uses MultiNest or PolyChord for nested sampling with C++/MPI parallelization and optional GPU acceleration.

## Quick Start

### Build
```bash
./autogen.sh && ./configure
make temponest && make temponest-install
make test
```

### Usage
```bash
tempo2 -gr temponest -f par.file tim.file -Cfile config.json
```

### Dependencies
- `$TEMPO2` environment variable
- MultiNest, GSL, BLAS/LAPACK, Eigen3, MPI
- Optional: ArrayFire for GPU acceleration

## Architecture

**Compositional Model System:**
```
model_t → model_space_t → [timing_model_t, efac_t, equad_t, ecorr_t, power_law_*_noise_t]
       → likelihood_t → temponest_v1_t
```

## Configuration

JSON-based with three main sections:
- `globals`: Analysis settings
- `sampler`: MultiNest parameters  
- `elements`: Model components

### Model Elements
- **Timing Model**: Tempo2 parameter handling (marginalise/fit modes)
- **EFAC/EQUAD**: Error scaling and white noise
- **ECORR**: Epoch correlations with multi-flag support
- **Power Law Noise**: Red/DM stochastic processes
- **Solar Wind**: Deterministic/stochastic electron density effects

## Tools

- **GUI**: `TempoNest_JSON.py` - Streamlit configuration interface
- **Plotting**: `plot_corner.py` - Enhanced corner plots with comparison mode

## Key Features

- **GPU Acceleration**: ArrayFire-based with 10-50x speedup potential
- **Robust Parameter Handling**: Automatic JUMP marginalization
- **Multi-flag ECORR**: Legacy-compatible epoch correlation modeling
- **Enhanced Solar Wind Support**: Both deterministic and stochastic models
- **Comparison Tools**: Multi-dataset posterior comparison capabilities

## Performance

- **MultiNest Live Points**: 500 (testing), 1000 (standard), 4000 (publication)
- **GPU Scaling**: Effective for >20 parameter problems
- **Resource Efficiency**: 1 GPU ≈ 300-400 CPU cores for high-dimensional problems

See `temponest_doc.md` for comprehensive documentation, examples, and technical details.