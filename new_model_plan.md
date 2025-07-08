# TempoNest New Chromatic Models Implementation Plan

## Overview
This document outlines the implementation plan for adding new chromatic signal models to TempoNest, based on Enterprise and Enterprise Extensions implementations. These models will enhance TempoNest's capability to handle various chromatic effects in pulsar timing data.

## Model Descriptions and Enterprise References

### 1. Chromatic Gaussian Signal

**Purpose**: Models transient chromatic events with Gaussian temporal profile (solar wind disturbances, ionospheric effects, ISM variations).

**Enterprise Reference**: 
- File: `reference_code/models_to_implement.notes` lines 21-38
- Function: `chrom_gaussian(toas, freqs, log10_Amp=-2.5, sign_param=1.0, t0=53890, sigma=81, idx=2)`

**Mathematical Form**:
```
Δt(t,ν) = sign × 10^(log10_Amp) × exp(-(t-t0)²/(2σ²)) × (1400/ν)^idx
```

**Enterprise Logic**:
```python
wf = 10**log10_Amp * np.exp(-(toas - t0)**2/2/sigma**2)
return np.sign(sign_param) * wf * (1400 / freqs) ** idx
```

**Parameters**:
- `log10_Amp`: Log amplitude (-10 to -5) [Enterprise: -10 to -5]
- `sign_param`: Sign of effect (-1 to 1) [Enterprise: -1 to 1]
- `t0`: Center time (MJD) [Enterprise: toas.min()/86400 to toas.max()/86400]
- `sigma`: Width in days [Enterprise: 0 to Tspan/2/86400]
- `idx`: Chromatic index (0-14) [Enterprise: 0 to 14]

**Enterprise Output Parameter Names**:
*Note: Enterprise typically uses GP representations rather than direct Gaussian signals*
- `<pulsar_name>_chromgauss_log10_Amp_chrom_gauss` - Log amplitude
- `<pulsar_name>_chromgauss_sign_param` - Sign parameter  
- `<pulsar_name>_chromgauss_t0` - Center time [MJD]
- `<pulsar_name>_chromgauss_sigma` - Width [days]
- `<pulsar_name>_chromgauss_idx` - Chromatic index (if variable)

**TempoNest Implementation Strategy**:
- Class: `chromatic_gaussian_t`
- Method: `apply()` modifies residuals vector directly
- Integration: Called during likelihood calculation before covariance operations

### 2. Chromatic Exponential Decay

**Purpose**: Models transient events with exponential decay (solar wind shocks, magnetospheric disturbances).

**Enterprise Reference**:
- File: `reference_code/enterprise_extensions/enterprise_extensions/chromatic/chromatic.py` lines 24-44
- Function: `chrom_exp_decay(toas, freqs, log10_Amp=-7, sign_param=-1.0, t0=54000, log10_tau=1.7, idx=2)`

**Mathematical Form**:
```
Δt(t,ν) = sign × 10^(log10_Amp) × H(t-t0) × exp(-(t-t0)/τ) × (1400/ν)^idx
```

**Enterprise Logic**:
```python
t0 *= const.day
tau = 10**log10_tau * const.day
ind = np.where(toas > t0)[0]
wf = 10**log10_Amp * np.heaviside(toas - t0, 1)
wf[ind] *= np.exp(-(toas[ind] - t0) / tau)
return np.sign(sign_param) * wf * (1400 / freqs) ** idx
```

**Parameters**:
- `log10_Amp`: Log amplitude (-10 to -5) [Enterprise: -10 to -2]
- `sign_param`: Sign of effect (-1 to 1) [Enterprise: -1 to 1]
- `t0`: Start time (MJD) [Enterprise: tmin to tmax]
- `log10_tau`: Log decay time in days (0 to 2.5) [Enterprise: 0 to 2.5]
- `idx`: Chromatic index (0-14) [Enterprise: 2 or vary 1-6]

**Enterprise Output Parameter Names**:
- `<pulsar_name>_<signal_name>_t0_dmexp` - Time of exponential minimum [MJD]
- `<pulsar_name>_<signal_name>_log10_Amp_dmexp` - Log10 amplitude of dip
- `<pulsar_name>_<signal_name>_log10_tau_dmexp` - Log10 of 1/e decay time [log10(days)]
- `<pulsar_name>_<signal_name>_sign_param` - Sign parameter (if varying)
- `<pulsar_name>_<signal_name>_idx` - Chromatic index (if variable)

**TempoNest Implementation Strategy**:
- Class: `chromatic_exponential_t`
- Method: `apply()` with Heaviside function and exponential decay
- Time conversion: Enterprise uses seconds, TempoNest uses MJD

### 3. DM/Chromatic Annual and Orbital Signals

**Purpose**: Models periodic variations due to Earth's orbital motion and binary orbital effects.

**Enterprise Reference**:
- File: `reference_code/models_to_implement.notes` lines 67-95
- Functions: `chrom_yearly_sinusoid()` and `chrom_orb_sinusoid()`

**Annual Signal**:
```python
def chrom_yearly_sinusoid(toas, freqs, log10_Amp=-7, phase=0, idx=2):
    wf = 10**log10_Amp * np.sin(2 * np.pi * const.fyr * toas + phase)
    return wf * (1400 / freqs) ** idx
```

**Orbital Signal**:
```python
def chrom_orb_sinusoid(toas, freqs, log10_Amp=-7, phase=0, idx=2, pb=365.25):
    forb = 1/(pb * 86400)
    wf = 10**log10_Amp * np.sin(2 * np.pi * forb * toas + phase)
    return wf * (1400 / freqs) ** idx
```

**Parameters**:
- `log10_Amp`: Log amplitude (-10 to -5) [Enterprise: -10 to -5]
- `phase`: Phase offset (0 to 2π) [Enterprise: 0 to 2π]
- `idx`: Chromatic index (typically 2) [Enterprise: 2]
- `pb`: Binary period in days (for orbital) [Enterprise: from pulsar parameters]

**Enterprise Output Parameter Names**:
- `<pulsar_name>_dm_s1yr_log10_Amp_dm1yr` - Log10 amplitude of annual signal
- `<pulsar_name>_dm_s1yr_phase_dm1yr` - Phase of annual sinusoid [0, 2π]
- `<pulsar_name>_dm_s1yr_idx` - Chromatic index (if variable)
- `<pulsar_name>_dmorb_log10_Amp_dmorb` - Log10 amplitude of orbital signal (if applicable)
- `<pulsar_name>_dmorb_phase_dmorb` - Phase of orbital sinusoid (if applicable)
- `<pulsar_name>_dmorb_pb` - Binary period [days] (if fitted)

**TempoNest Implementation Strategy**:
- Class: `chromatic_periodic_t`
- Support both annual and orbital modes
- Annual frequency: 1/(365.25*86400) Hz
- Orbital frequency: 1/(pb*86400) Hz where pb from pulsar parameters

### 4. Chromatic GP Noise Model

**Purpose**: Models general stochastic chromatic noise with power-law spectrum and frequency dependence.

**Enterprise Reference**:
- File: `reference_code/enterprise_extensions/enterprise_extensions/blocks.py` lines 863-1166
- Function: `chromatic_noise_block()` with `gp_kernel="diag"` and `psd="powerlaw"`

**Mathematical Form**:
```
Power Spectral Density: P(f) = A² × f^(-γ)
Frequency scaling: ∝ (ν_ref/ν)^idx
```

**Enterprise Logic**:
```python
# Diagonal chromatic GP with power-law PSD
chrom_gp = gp_signals.FourierBasisGP(
    powerlaw(log10_A=log10_A_chrom, gamma=gamma_chrom),
    components=components,
    Tspan=Tspan,
    name='chrom_gp'
)
```

**Parameters**:
- `log10_A`: Log amplitude (-20 to -11) [Enterprise: -20 to -11]
- `gamma`: Spectral index (0 to 7) [Enterprise: 0 to 7]
- `idx`: Chromatic frequency index (0 to 14) [Enterprise: 0 to 14, default 4]

**Enterprise Output Parameter Names**:
- `<pulsar_name>_chrom_gp_log10_A` - Log10 amplitude of chromatic GP
- `<pulsar_name>_chrom_gp_gamma` - Spectral index of chromatic GP
- `<pulsar_name>_chrom_gp_idx` - Chromatic frequency index
- `<pulsar_name>_chrom_gp_<mode>` - Individual Fourier coefficients (mode = 0,1,2...)

**TempoNest Implementation Strategy**:
- Class: `chromatic_gp_noise_t` (similar to existing `pl_red_noise_t`)
- Method: Add chromatic frequency scaling to existing GP framework
- Integration: Modify covariance matrix with frequency-dependent scaling
- Basis: Use same Fourier basis as red noise but with chromatic scaling

### 5. Enhanced GP Solar Wind Model

**Purpose**: Sophisticated solar wind modeling using Gaussian Process with Fourier basis.

**Enterprise Reference**:
- File: `reference_code/enterprise_extensions/enterprise_extensions/chromatic/solar_wind.py` lines 146-184
- Function: `createfourierdesignmatrix_solar_dm()`

**Enterprise Logic**:
```python
# Get base fourier design matrix
F, Ffreqs = utils.createfourierdesignmatrix_red(toas, nmodes=nmodes, ...)
# Calculate solar wind geometric factors
theta, R_earth, _, _ = theta_impact(planetssb, sunssb, pos_t)
dm_sol_wind = dm_solar(1.0, theta, R_earth)
dt_DM = dm_sol_wind * 4.148808e3 / (freqs**2)
# Scale Fourier basis by solar wind geometry and frequency
return F * dt_DM[:, None], Ffreqs
```

**Enterprise GP Setup**:
```python
log10_A_sw = parameter.Uniform(-10, -2)  
gamma_sw = parameter.Uniform(0, 4)
sw_prior = utils.powerlaw(log10_A=log10_A_sw, gamma=gamma_sw)
sw_components = 30
sw_basis = createfourierdesignmatrix_solar_dm(nmodes=sw_components, Tspan=Tspan)
sw = gp_signals.BasisGP(sw_prior, sw_basis, name='gp_sw')
```

**Enterprise Output Parameter Names**:

**Deterministic Solar Wind:**
- `<pulsar_name>_n_earth_n_earth` - Solar wind electron density at 1 AU [cm^-3]

**GP Solar Wind (Power-law basis):**
- `<pulsar_name>_gp_sw_log10_A_sw` - Log10 amplitude of solar wind GP
- `<pulsar_name>_gp_sw_gamma_sw` - Spectral index of solar wind GP
- `<pulsar_name>_gp_sw_<mode>` - Individual Fourier coefficients (mode = 0,1,2...)

**GP Solar Wind (Periodic kernel):**
- `<pulsar_name>_gp_sw_log10_sigma` - Log10 variance amplitude
- `<pulsar_name>_gp_sw_log10_ell` - Log10 correlation length scale [log10(days)]
- `<pulsar_name>_gp_sw_log10_p` - Log10 periodicity scale [log10(years)]
- `<pulsar_name>_gp_sw_log10_gam_p` - Log10 periodicity strength

**GP Solar Wind (Squared-exponential kernel):**
- `<pulsar_name>_gp_sw_log10_sigma` - Log10 variance amplitude
- `<pulsar_name>_gp_sw_log10_ell` - Log10 correlation length scale [log10(days)]

**TempoNest Implementation Strategy**:
- Enhance existing `stochastic_solar_wind_t`
- Add Fourier basis matrix calculation
- Implement power-law spectral prior
- Add geometric solar wind scaling factors
- Integrate with existing noise covariance framework

## Implementation Architecture

### Current TempoNest Model Element Structure

**Base Class**: `model_element_t`
- Parameter management via `parameter_map_`
- Virtual methods: `is_valid_parameter()`, `is_fully_specified()`, `get_name()`, `set_parameter()`, `write_to_par_file()`

**Existing Examples**:
- `deterministic_solar_wind_t`: Simple linear correction to residuals
- `stochastic_solar_wind_t`: Adds noise variance terms
- `pl_red_noise_t`: Power-law red noise with Fourier frequencies

### Integration Points

1. **Likelihood Calculation** (`likelihood_t::operator()`):
   - Deterministic models: Apply corrections to residuals before covariance operations
   - Stochastic models: Modify noise covariance matrix

2. **Parameter Management** (`model_space_t`):
   - Each element manages its own parameters
   - Parameters indexed for sampling algorithms

3. **JSON Configuration**:
   - Model elements configured via JSON files
   - Parameter priors and ranges specified

## Implementation Plan

### Phase 1: Deterministic Chromatic Models (Week 1)

#### 1.1 Chromatic Gaussian Signal
- **File**: `src/core/models/elements/chromatic_gaussian.h/.cpp`
- **Parameters**: log10_Amp, sign_param, t0, sigma, idx
- **Apply Method**: Gaussian waveform with frequency scaling
- **Validation**: Compare against Enterprise `chrom_gaussian` output

#### 1.2 Chromatic Exponential Decay
- **File**: `src/core/models/elements/chromatic_exponential.h/.cpp`
- **Parameters**: log10_Amp, sign_param, t0, log10_tau, idx
- **Apply Method**: Heaviside × exponential with frequency scaling
- **Validation**: Compare against Enterprise `chrom_exp_decay` output

#### 1.3 DM Periodic Signals
- **File**: `src/core/models/elements/chromatic_periodic.h/.cpp`
- **Parameters**: log10_Amp, phase, idx, pb (optional)
- **Modes**: Annual (fixed frequency) and orbital (pulsar-dependent frequency)
- **Apply Method**: Sinusoidal waveform with frequency scaling
- **Validation**: Compare against Enterprise periodic functions

#### 1.4 Chromatic GP Noise
- **File**: `src/core/models/elements/chromatic_gp_noise.h/.cpp`
- **Parameters**: log10_A, gamma, idx
- **Method**: Enhance existing GP framework with chromatic frequency scaling
- **Integration**: Modify covariance matrix calculation
- **Validation**: Compare against Enterprise chromatic_noise_block output

### Phase 2: Dynamic Scaling Implementation (Hybrid Approach)

#### 2.1 Mathematical Foundation

**Chromatic GP Model:**
The chromatic Gaussian Process has separable temporal and frequency components:
```
Chromatic GP = F_temporal(t) ⊗ F_frequency(ν)
```

Where:
- **Temporal Component**: `F_temporal = [cos(2πf₁t), sin(2πf₁t), cos(2πf₂t), sin(2πf₂t), ...]`
- **Frequency Component**: `F_frequency = (ν_ref/ν)^idx` where `ν_ref = 1400 MHz`

**Power Spectral Density:**
```
P(f_temporal) = A² × f_temporal^(-γ)
```

**Full Covariance Matrix:**
```
C_chromatic = F_chromatic × P × F_chromatic^T
```

#### 2.2 Dynamic Scaling Challenge

**Problem**: When `idx` is a fitted parameter, the frequency scaling `(ν_ref/ν)^idx` changes on every likelihood evaluation, but the design matrix `F_chromatic` is constructed once at startup.

**Traditional Solutions:**
1. **Matrix Recalculation**: Update design matrix every likelihood call - `O(nobs × nfreqs)` operations
2. **Deferred Scaling**: Apply scaling in covariance computation - requires major architecture changes

#### 2.3 Hybrid Approach Solution

**Key Insight**: The effect of changing `idx` is mathematically equivalent to scaling the power coefficients by the ratio of integrated frequency scaling factors.

**Mathematical Derivation:**

If the design matrix was constructed with reference index `idx_ref`:
```
F_ref = F_temporal ⊗ (ν_ref/ν)^idx_ref
```

But we want the effect of current index `idx_curr`:
```
F_curr = F_temporal ⊗ (ν_ref/ν)^idx_curr
```

The ratio of squared scaling factors (which affects power) is:
```
correction_factor = Σₖ[(ν_ref/νₖ)^(2×idx_curr)] / Σₖ[(ν_ref/νₖ)^(2×idx_ref)]
```

**Implementation Steps:**

1. **Design Matrix Construction** (once at startup):
   ```cpp
   double reference_idx = 4.0;  // Fixed reference for matrix construction
   for (int k = 0; k < nobs; k++) {
       ChromVec[k] = pow(ref_freq / obs_freq[k], reference_idx);
       // Build temporal Fourier basis scaled by reference chromatic index
   }
   ```

2. **Dynamic Scaling Correction** (every likelihood call):
   ```cpp
   double calculate_scaling_correction(double current_idx, double reference_idx) {
       double sum_current = 0.0, sum_reference = 0.0;
       for (int k = 0; k < nobs; k++) {
           double freq_ratio = ref_freq / obs_freq[k];
           sum_current += pow(freq_ratio, 2 * current_idx);
           sum_reference += pow(freq_ratio, 2 * reference_idx);
       }
       return sum_current / sum_reference;
   }
   ```

3. **Power Coefficient Update**:
   ```cpp
   if (is_idx_fitted()) {
       double current_idx = get_parameter_value("idx");
       double correction_factor = calculate_scaling_correction(current_idx, 4.0);
       power_coefficients *= correction_factor;
   }
   ```

#### 2.4 Scientific Validation Strategy

**Test 1: DM-Noise-Only Convergence Test**
- **Dataset**: Simulated data with only DM variations (idx should converge to ~2)
- **Configuration**: Fit chromatic GP with `idx` uniform prior [0, 3]
- **Expected Result**: Posterior should peak around idx ≈ 2.0
- **Validation**: Confirms proper frequency scaling behavior

**Test 2: Enterprise Comparison**
- **Method**: Compare TempoNest vs Enterprise for identical configurations
- **Parameters**: Fixed `idx` values (2, 4, 6) and variable `idx`
- **Metrics**: Log-likelihood values, parameter posteriors
- **Tolerance**: `|L_TempoNest - L_Enterprise| < 1e-10`

**Test 3: Fixed vs Variable idx Consistency**
- **Method**: Run with `idx` fixed at specific value vs variable with tight prior around same value
- **Expected Result**: Should produce nearly identical results
- **Validation**: Ensures dynamic scaling implementation correctness

#### 2.5 Performance Benefits

**Computational Complexity:**
- **Matrix Recalculation**: `O(nobs × nfreqs)` ≈ 30,000 operations per likelihood
- **Hybrid Approach**: `O(nobs)` ≈ 1,000 operations per likelihood
- **Speedup**: ~30× faster per likelihood evaluation

**Memory Efficiency:**
- No additional matrix storage required
- Minimal memory overhead for scaling correction
- Maintains current TempoNest architecture

#### 2.6 Implementation Details

**Code Modifications Required:**

1. **chromatic_gp_noise.cpp**: Add `calculate_scaling_correction()` method
2. **chromatic_gp_noise.cpp**: Enhance `apply()` method with dynamic scaling
3. **model_space.cpp**: Update initialization message for variable idx
4. **JSON configuration**: Support for variable idx mode

**Backward Compatibility:**
- Fixed `idx` mode unchanged (no performance impact)
- Variable `idx` mode uses hybrid approach
- All existing configurations continue to work

### Phase 3: Enhanced GP Solar Wind (Week 3)

#### 3.1 Fourier Basis Implementation
- Enhance `stochastic_solar_wind_t` with Fourier basis matrix
- Add `calculate_fourier_basis()` method
- Implement solar wind geometric scaling factors

#### 3.2 Power-Law Prior Integration
- Add power-law spectral density calculation
- Integrate with existing covariance matrix framework
- Add `log10_A` and `gamma` parameters

#### 3.3 Solar Wind Geometry
- Implement `theta_impact()` equivalent for solar angle calculation
- Add Earth position and pulsar position handling
- Scale Fourier basis by geometric factors

### Phase 4: Model Integration (Week 4)

#### 4.1 Update Core Headers
- Add new model includes to `model_element.h`
- Update model factory/instantiation code
- Ensure proper parameter indexing

#### 4.2 JSON Configuration Support
- Add JSON parsing for new model parameters
- Create example configuration files
- Update parameter validation

#### 4.3 Likelihood Integration
- Ensure new models called in proper order
- Validate deterministic vs stochastic model separation
- Test parameter gradients for optimization

### Phase 5: Testing and Validation (Week 5)

#### 5.1 Unit Tests
- Create test cases for each model
- Compare outputs with Enterprise implementations
- Validate parameter ranges and priors

#### 5.2 Integration Tests
- Test full likelihood calculations
- Validate sampling convergence
- Performance benchmarking

#### 5.3 Example Configurations
- Create realistic test cases
- Document parameter selection guidelines
- Provide usage examples

## Technical Implementation Details

### Time Unit Conversions
- **Enterprise**: Uses seconds for TOAs, converts MJD×86400
- **TempoNest**: Uses MJD directly
- **Conversion Strategy**: Convert internally to match Enterprise formulas

### Frequency Scaling
- **Standard**: `(1400/freq)^idx` where 1400 MHz is reference
- **DM-like**: idx=2
- **Scattering-like**: idx=4
- **Variable**: idx can be fitted parameter

### Parameter Prior Ranges
Follow Enterprise conventions:
- `log10_Amp`: Uniform(-10, -5) for most signals
- `phase`: Uniform(0, 2π)
- `t0`: Uniform(tmin, tmax) based on data span
- `sigma`: Uniform(0, Tspan/2) in days
- `log10_tau`: Uniform(0, 2.5) in log10(days)
- `idx`: Uniform(0, 14) when variable

### Code Style and Conventions
- Follow existing TempoNest C++ style
- Use Eigen for vector operations
- Maintain const-correctness
- Include comprehensive error checking
- Add detailed comments referencing Enterprise equivalents

## Validation Strategy

1. **Mathematical Verification**: Compare analytical formulas
2. **Numerical Validation**: Test identical inputs produce identical outputs
3. **Parameter Coverage**: Test full prior ranges
4. **Edge Cases**: Test boundary conditions and special cases
5. **Integration Testing**: Validate within full likelihood framework

## Expected Outcomes

After implementation, TempoNest will have:
1. Enhanced chromatic signal modeling capability
2. Better compatibility with Enterprise analysis workflows
3. More sophisticated solar wind modeling
4. Improved capability for detecting and characterizing transient events
5. Foundation for future GP-based signal implementations

This implementation will significantly expand TempoNest's scientific capabilities while maintaining its computational efficiency advantages.

## Enterprise Parameter Naming Conventions

Enterprise follows a systematic naming convention for all output parameters. Understanding this is crucial for compatibility and validation.

### General Pattern
```
<pulsar_name>_<signal_name>_<parameter_name>
```

**Components:**
- `<pulsar_name>`: Pulsar identifier (e.g., "J0437-4715", "J1853+1303")
- `<signal_name>`: Signal identifier (e.g., "dm_s1yr", "gp_sw", "dmexp")
- `<parameter_name>`: Specific parameter (e.g., "log10_Amp", "gamma", "t0")

### Parameter Type Conventions

**Amplitude Parameters:**
- Always in log10 space: `log10_A`, `log10_Amp`, `log10_sigma`
- Typical ranges: -10 to -2 for most signals

**Time Parameters:**
- Times in MJD: `t0`, `t_init`, `t_final`
- Time scales in log10 days: `log10_tau`, `log10_ell`
- Periods in days: `pb` (binary period)

**Spectral Parameters:**
- Power-law indices: `gamma` (typically 0-7)
- Chromatic indices: `idx` (typically 0-14)

**Phase Parameters:**
- Always in radians: `phase` (0 to 2π)

**GP Kernel Parameters:**
- `log10_sigma`: Log10 variance amplitude
- `log10_ell`: Log10 correlation length scale
- `log10_p`: Log10 periodicity scale
- `log10_gam_p`: Log10 periodicity strength

### Signal Name Conventions

**Common Enterprise Signal Names:**
- `dm_s1yr`: DM annual signal
- `dmorb`: DM orbital signal  
- `dmexp`: DM exponential dip
- `dm_cusp`: DM exponential cusp
- `gp_sw`: Solar wind Gaussian process
- `n_earth`: Deterministic solar wind
- `red_noise`: Red noise process
- `dm_gp`: DM Gaussian process
- `chrom_gp`: Chromatic Gaussian process

### TempoNest Parameter Mapping Strategy

To ensure compatibility, TempoNest parameter names should follow Enterprise conventions:

**Recommended TempoNest Parameter IDs:**
```cpp
// Chromatic Gaussian
"log10_Amp_chrom_gauss"
"sign_param" 
"t0"
"sigma"
"idx"

// Chromatic Exponential  
"log10_Amp_dmexp"
"sign_param"
"t0_dmexp" 
"log10_tau_dmexp"
"idx"

// DM Annual/Orbital
"log10_Amp_dm1yr"
"phase_dm1yr"
"idx"
"pb" (for orbital)

// Enhanced Solar Wind GP
"log10_A_sw" 
"gamma_sw"
"log10_sigma_sw"
"log10_ell_sw"
```

### Output Compatibility

TempoNest output files should use parameter names that can be directly compared with Enterprise results:

**PAR File Output:**
```
# Use Enterprise-compatible parameter names in comments
# Chromatic Gaussian: J1234+5678_chromgauss_log10_Amp_chrom_gauss = -6.234 +/- 0.123
```

**Chain Files:**
```
# Column headers should match Enterprise naming when possible
# log10_Amp_chrom_gauss  sign_param  t0  sigma  idx
```

This naming consistency will enable:
1. Direct comparison with Enterprise results
2. Easy validation of implementations  
3. Smooth workflow transitions between codes
4. Consistent scientific interpretation

## Target Parameter Coverage Analysis

Based on the J1909-3744 parameter list, here's our implementation coverage:

### ✅ ALREADY IMPLEMENTED in TempoNest:
1. `J1909-3744_KAT_MKBF_efac` → **EFAC** (`efac_t` class)
2. `J1909-3744_KAT_MKBF_log10_ecorr` → **ECORR** (`ecorr_t` class)  
3. `J1909-3744_KAT_MKBF_log10_tnequad` → **EQUAD** (`equad_t` class)
4. `J1909-3744_red_noise_gamma` → **Red Noise** (`pl_red_noise_t` class)
5. `J1909-3744_red_noise_log10_A` → **Red Noise** (`pl_red_noise_t` class)
6. `J1909-3744_dm_gp_gamma` → **DM Noise** (`pl_dm_noise_t` class - spectral_index parameter)
7. `J1909-3744_dm_gp_log10_A` → **DM Noise** (`pl_dm_noise_t` class - amplitude parameter)

### ✅ PLANNED FOR IMPLEMENTATION (This Project):
8. `J1909-3744_chrom_gp_gamma` → **Chromatic GP Noise** (`chromatic_gp_noise_t`)
9. `J1909-3744_chrom_gp_idx` → **Chromatic GP Noise** (`chromatic_gp_noise_t`)
10. `J1909-3744_chrom_gp_log10_A` → **Chromatic GP Noise** (`chromatic_gp_noise_t`)
11. `J1909-3744_sw_gp_gamma` → **Enhanced Solar Wind GP** (enhance `stochastic_solar_wind_t`)
12. `J1909-3744_sw_gp_log10_A` → **Enhanced Solar Wind GP** (enhance `stochastic_solar_wind_t`)

### ❌ NOT NEEDED:
13. `nmodel` → **Model Selection Parameter** (excluded per user request)

## Coverage Assessment: **92% Complete** (12/13 parameters)

Excellent news! TempoNest already has much better coverage than initially assessed:

### Current Status:
- **Already implemented**: 7/13 parameters (54%)
- **Will be implemented**: 5/13 parameters (38%) 
- **Not needed**: 1/13 parameters (8%)

### **Total coverage after this project: 12/12 needed parameters = 100%!**

## Final Recommendation:

The current implementation plan (Phases 1-5) will achieve **complete coverage** of all the scientifically relevant parameters from your J1909-3744 analysis. TempoNest already has excellent noise modeling capabilities, and adding the chromatic models will provide full compatibility with Enterprise analyses.

**Key insight**: TempoNest's existing `pl_dm_noise_t` already handles DM GP noise, so the focus should be on the new chromatic signal capabilities that aren't currently available.