# TempoNest Debugging Guide

This guide explains how to use TempoNest's verbosity settings to debug issues and monitor sampling progress.

## Verbosity Levels

TempoNest supports three verbosity levels to balance between information and performance:

1. **Quiet** (`verbose: false`): Minimal output - only critical errors and final results
2. **Normal** (`verbose: true`): Progress updates every 100 iterations and important milestones  
3. **Full** (`verbose: "full"`): Complete debug output including all computational details

## Setting Verbosity

Set the verbosity level in your sampler configuration file (e.g., `sampler_polychord.json`):

```json
{
    "sampler": "polychord",
    "verbose": true,    // Options: false, true, "full"
    "num_live": 200,
    ...
}
```

### Examples:

**Quiet mode** - for production runs where you only need results:
```json
"verbose": false
```
Output: Only critical errors and completion message

**Normal mode** - recommended for monitoring progress:
```json
"verbose": true
```
Output:
- "PolyChord sampling started - first dead point processed"
- Progress updates every 100 iterations
- Final results

**Full debug mode** - for troubleshooting:
```json
"verbose": "full"
```
Output: All Normal mode output plus detailed debug information (see below)

## Understanding the Debug Output

When verbose mode is enabled, TempoNest provides detailed step-by-step information about the likelihood calculation process. Here's what each debug step means:

### Parameter Information

```
DEBUG STEP 0 - Parameter values [0:10]: 
  param[0] = -14.714957
  param[1] = 4.636625
  ...
  ... (total 76 parameters)
```

**Interpretation**: Shows the first 10 parameter values and total parameter count. These are the transformed parameters from the sampler (PolyChord/MultiNest) in physical units.

### Timing Residuals (Step 1)

```
DEBUG STEP 1 - Residuals: min=-0.000046, max=0.000030, has_nan=NO
```

**Interpretation**: 
- Shows the range of timing residuals after applying timing model parameters
- `has_nan=NO` confirms no numerical issues in residual calculation
- Values are in seconds

### White Noise Components (Steps 2-4)

```
DEBUG STEP 2 - Initial noise: min=0.000000, max=0.000031, has_nan=NO
DEBUG STEP 3 - After EFAC: min=0.000000, max=0.000000, has_nan=NO
DEBUG STEP 4 - After EQUAD: min=0.000000, max=0.000001, has_nan=NO
```

**Interpretation**:
- **Step 2**: Initial TOA uncertainties from the data files
- **Step 3**: After applying EFAC (multiplicative white noise scaling)
- **Step 4**: After adding EQUAD (additive white noise in quadrature)
- Values represent noise variance (σ²)

### Solar Wind Effects (Step 5)

```
DEBUG STEP 5 - After Det Solar Wind: min=-0.000046, max=0.000030, has_nan=NO
```

**Interpretation**: Shows residuals after applying deterministic solar wind delays (if included in the model).

### Noise Floor (Step 6.5)

```
DEBUG STEP 6.5 - After noise floor: min=0.000000, max=0.000001, has_nan=NO
```

**Interpretation**: A minimum noise floor (1e-12 seconds) is applied to prevent division by zero in subsequent calculations.

### Noise Inverse (Step 7)

```
DEBUG STEP 7 - After noise inverse: min=2330840.075968, max=1000000000000.000000, has_nan=NO
```

**Interpretation**: Shows the inverted noise values (1/σ²). Large values indicate very precise measurements.

### Red Noise Coefficients (Step 7.1)

```
DEBUG STEP 7.1 - Red Noise applied: coeffs=60, freq_det_contrib=-2624.642277
DEBUG STEP 7.1 - Red powercoeff range: min=0.000000, max=0.000000
DEBUG STEP 7.1 - red_powercoeff[0:5]: 0.000000 0.000000 0.000000 0.000000 0.000000
```

**Interpretation**:
- `coeffs=60`: Number of Fourier coefficients for red noise
- `freq_det_contrib`: Contribution to the frequency domain determinant
- `powercoeff`: Power spectral coefficients (may show as 0.000000 due to display precision for very small values)

### DM Noise Coefficients (Step 7.2)

```
DEBUG STEP 7.2 - DM Noise applied: coeffs=200, freq_det_contrib=-8210.222948
DEBUG STEP 7.2 - DM powercoeff range: min=0.000000, max=0.000000
```

**Interpretation**: Similar to red noise but for dispersion measure (DM) variations.

### ECORR Coefficients (Step 7.3)

```
DEBUG STEP 7.3 - ECORR applied: coeffs=3798, freq_det_contrib=-114253.810323
DEBUG STEP 7.3 - ECORR powercoeff range: min=0.000000, max=0.000000
```

**Interpretation**:
- `coeffs=3798`: Number of ECORR epochs across all backends
- Large `freq_det_contrib` is normal due to many epochs

### Power Coefficient Summary (Step 7.4)

```
DEBUG STEP 7.4 - POWERCOEFF SUMMARY: total_coeffs=4058, overall_min=0.000000, overall_max=0.000000
DEBUG STEP 7.4 - Extreme values (>1e15): count=0
```

**Interpretation**:
- Total number of noise coefficients (red + DM + ECORR)
- Checks for numerical instabilities (extreme values)
- Values shown as 0.000000 may actually be very small (e.g., 1e-16)

### Likelihood Components (Step 8)

```
DEBUG STEP 8 - Time domain: timelike=17388.247034 (finite=YES), tdet=-146355.038614 (finite=YES)
DEBUG STEP 8.0 - Using GPU path, optimized=YES
DEBUG STEP 8.4 - GPU result: likelihood=65307.270061 (finite=YES)
```

**Interpretation**:
- `timelike`: Time-domain chi-squared contribution
- `tdet`: Time-domain determinant (log of noise terms)
- Shows whether GPU acceleration is being used
- Final likelihood value should be finite

### Final Result (Step 9)

```
DEBUG STEP 9 - Final likelihood: 65307.270061 (finite=YES)
DEBUG COMPONENTS - tdet=-146355.038614, freq_det=-114253.810323, timelike=17388.247034, uniform_prior=0.000000
```

**Interpretation**: 
- Final log-likelihood value and its components
- Should always be finite for valid parameter combinations

## Common Issues and What to Look For

### 1. NaN Values
If any step shows `has_nan=YES`, this indicates numerical issues:
- Check parameter bounds in configuration
- Look for extreme parameter values in Step 0
- Verify data file integrity

### 2. Extreme Coefficients
If Step 7.4 shows `Extreme values (>1e15): count > 0`:
- May indicate parameter bounds are too wide
- Check for very small noise parameters causing numerical overflow

### 3. Non-finite Likelihood
If Step 9 shows `(finite=NO)`:
- Usually caused by invalid parameter combinations
- Check for division by zero or log of negative numbers
- Review parameter priors and bounds

### 4. GPU vs CPU Differences
Step 8.0 shows computation path:
- `Using GPU path, optimized=YES`: Most efficient
- `Using GPU path, optimized=NO`: Standard GPU
- `Using CPU path`: Fallback when GPU unavailable

## Parameter Counting

The debug output helps verify correct parameter counting:

```
PARAM COUNT: Power Law Red Noise added 2 parameters, param_index: 0 -> 2
PARAM COUNT: Power Law DM Noise added 2 parameters, param_index: 2 -> 4
PARAM COUNT: EFAC added 20 parameters, param_index: 4 -> 24
PARAM COUNT: EQUAD added 20 parameters, param_index: 24 -> 44
PARAM COUNT: ECORR added 31 parameters, param_index: 44 -> 75
PARAM COUNT: Deterministic Solar Wind added 1 parameters, param_index: 75 -> 76
```

This shows the sequential assignment of parameter indices and helps identify indexing issues.

## Performance Considerations

Verbose mode adds overhead to calculations. For production runs:
- Disable verbose mode once debugging is complete
- Use verbose mode only for initial test runs
- Consider enabling only for specific problematic pulsars

## Integration with Samplers

When using PolyChord or MultiNest, verbose output appears interspersed with sampler output:
```
PolyChord prior transform #0: cube→physical: 0.321261→-14.714957 ...
PolyChord likelihood eval #0: params[0:4]=-14.714957 4.636625 ...
DEBUG STEP 0 - Parameter values [0:10]: ...
...
PolyChord likelihood result #0: loglike=59895.182613
```

This helps correlate specific parameter combinations with likelihood calculation issues.

## Tips for Effective Debugging

1. **Start with a small test**: Use verbose mode on a single pulsar first
2. **Check parameter ranges**: Ensure all parameters have reasonable bounds
3. **Monitor for patterns**: Look for systematic issues across likelihood evaluations
4. **Compare with known good configurations**: Use reference implementations to validate
5. **Disable GPU initially**: Test with CPU first to isolate GPU-specific issues

## Related Documentation

- [Troubleshooting Guide](troubleshooting.md) - Common issues and solutions
- [Configuration Guide](configuration.md) - Parameter settings
- [GPU Optimization](gpu_optimization.md) - GPU-specific debugging