# TempoNest GPU Optimization and PolyChord Integration Progress

## Previous Work Completed
- Successfully implemented GPU optimization using ArrayFire
- Achieved 30% performance improvement with GPU acceleration
- GPU optimization is functional and providing expected speedup

## Current Issue: PolyChord Integration Problem

### Problem Statement
PolyChord sampler hangs after printing "started sampling" message. No sampling progress occurs, and no chain files are generated. This is critical because **most datasets are >30 dimensions, which MultiNest is not good for**.

### Investigation Summary

#### Initial Hypothesis: Parameter Configuration
- **Attempted Fix 1**: Reduced sampling parameters
  - `num_live: 500 → 100`, `precision_criterion: 0.5 → 0.1`
  - **Result**: Still hung at "started sampling"

- **Attempted Fix 2**: Added debug logging
  - Confirmed likelihood wrapper is being called during live point generation
  - **Discovery**: Live points generated successfully, but main sampling loop never progresses

- **Attempted Fix 3**: Further parameter reduction
  - `num_live: 100 → 25`, `num_repeats: 1`
  - **Result**: Same hanging behavior

#### Root Cause Analysis: Parameter Validation Issues
- **Discovery**: `num_repeats` must be ≥ 5×nDims for slice sampling
- **Missing Parameters**: `logzero`, `nprior`, `nfail` not set (required by PolyChord)
- **Implemented Fix**: Automatic correction in `polychord.cpp`:
  ```cpp
  // Auto-correct num_repeats if too low
  int recommended_repeats = 5 * pc_settings.nDims;
  pc_settings.num_repeats = std::max(settings.num_repeats, recommended_repeats);
  
  // Set missing required parameters
  pc_settings.logzero = -1e30;
  pc_settings.nprior = -1;
  pc_settings.nfail = -1;
  ```
- **Result**: Parameters corrected but **still hangs at "started sampling"**

#### MPI Compilation Issue Investigation
- **Discovery**: PolyChord compiled with `MPI=1` by default
- **Theory**: MPI initialization conflict causing deadlock in single-threaded context
- **Attempted Fix**: Rebuilt PolyChord with `MPI=0`
  - Changed `external/PolyChordLite/Makefile` line 18: `MPI=1` → `MPI=0`
  - Cleaned and rebuilt PolyChord library without MPI support
  - Rebuilt and reinstalled TempoNest
- **Result**: **Still hangs at "started sampling"**

### Current Status
- ✅ GPU optimization working (30% speedup)
- ✅ Parameter validation fixes implemented
- ✅ MPI compilation issue resolved
- ❌ **PolyChord still hangs during main sampling loop**

### Technical Details Investigated

#### PolyChord Code Analysis
- **Hang Location**: After "started sampling" message in `feedback.f90`
- **Main Loop**: `nested_sampling.F90` enters while loop but never progresses
- **Live Points**: Successfully generated (25 points visible in `_phys_live.txt`)
- **Likelihood Calls**: Working correctly during initialization phase

#### Compilation Flags Verified
- Non-MPI build uses `gfortran`, `gcc`, `g++` (not `mpifort`, `mpicc`, `mpicxx`)
- Optimization flags: `-Ofast` for performance
- No MPI preprocessor directives active (`-DMPI` and `-DUSE_MPI` removed)

### Next Steps for Investigation
The fundamental issue remains unresolved. Potential areas to investigate:

1. **Threading Issues**: PolyChord may have OpenMP conflicts
2. **Library Dependencies**: Missing mathematical libraries (LAPACK, BLAS)
3. **Fortran Runtime**: Potential issues with Fortran-C++ interface
4. **PolyChord Version**: May need different PolyChord version or build configuration
5. **Algorithm Logic**: Possible infinite loop in slice sampling implementation

### Files Modified
- `external/PolyChordLite/Makefile`: `MPI=1` → `MPI=0`
- `src/core/samplers/polychord.cpp`: Added parameter validation and auto-correction
- `tests/data/tnest_config_test/sampler_polychord.json`: Various parameter configurations tested

### Configuration Used
```json
{
    "sampler": {
        "type": "polychord",
        "output_root": "results_polychord/TNest-",
        "sample": true,
        "num_live": 25,
        "num_repeats": 1,
        "precision_criterion": 0.5,
        "do_clustering": false,
        "feedback": 2,
        "posteriors": true,
        "equals": true,
        "write_resume": true,
        "write_paramnames": true,
        "read_resume": false,
        "write_live": true,
        "write_dead": true
    }
}
```

**Note**: `num_repeats=1` is automatically corrected to `5×nDims` by the validation code.

### Priority
**HIGH** - This is blocking the ability to analyze high-dimensional datasets (>30D) effectively.