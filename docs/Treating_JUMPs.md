# JUMP-Type Parameter Treatment in Tempo2 and TempoNest

## Overview

This document provides comprehensive documentation on how JUMP-type parameters (JUMP, FDJUMP, FDJUMPDM) are handled in Tempo2 versus TempoNest, covering the conceptual, mathematical, and implementation differences between the two systems.

## Parameter Types and Physical Meaning

### 1. JUMP Parameters
**Physical Purpose**: Correct for systematic timing offsets between different observing systems, telescopes, or data processing pipelines.

**Mathematical Model**:
```
δt_jump = constant_offset
```

**Typical Use Cases**:
- Different telescope systems with systematic timing differences
- Backend/receiver changes that introduce timing offsets
- Data processing pipeline differences

**Parameter File Format**:
```
JUMP -sys EFF.EBPP.1360 0.00059159454699566 1
JUMP -f 430_ASP 0.00032357032081991 1
JUMP -pta MPTA 0.0012619966023039 1
```

### 2. FDJUMP Parameters (Frequency-Dependent Jumps)
**Physical Purpose**: Model frequency-dependent timing offsets with configurable power-law dependence.

**Mathematical Model**:
```
δt_fdjump = amplitude × (f/f_ref)^index
```
Where:
- `amplitude`: Fitted jump amplitude
- `f`: Observation frequency
- `f_ref`: Reference frequency (typically 1 GHz)
- `index`: User-specified frequency dependence (1, 2, 3, ...)

**Frequency Dependencies**:
- **FDJUMP1**: `δt ∝ f^1` (linear frequency dependence)
- **FDJUMP2**: `δt ∝ f^2` (quadratic frequency dependence)  
- **FDJUMP3**: `δt ∝ f^3` (cubic frequency dependence)

**Parameter File Format**:
```
FDJUMP1 -pta PPTA 1.6225418269989e-05 1      # First-order (f^1)
FDJUMP2 -pta PPTA -1.9380154188785e-05 1     # Second-order (f^2)
FDJUMP3 -pta PPTA 8.413771226843e-06 1       # Third-order (f^3)
```

### 3. FDJUMPDM Parameters (DM-like Jumps)
**Physical Purpose**: Model DM-like (dispersion measure) systematic offsets between observing systems.

**Mathematical Model**:
```
δt_fdjumpdm = amplitude × (f/f_ref)^(-2) / DM_CONST
```
This mimics the frequency dependence of interstellar dispersion: `δt ∝ f^(-2)`.

**Physical Interpretation**:
- Systematic differences in DM measurement/correction between systems
- Instrumental dispersion effects
- Frequency-dependent calibration errors

**Parameter File Format**:
```
FDJUMPDM -pta NANOGrav 0.021768097764128 1   # DM-like (f^-2) jump
```

## Tempo2 Implementation

### Conceptual Approach
**Philosophy**: "Fit all user-specified parameters to minimize timing residuals"

Tempo2 treats all JUMP-type parameters as normal fitted parameters in a least-squares framework. Users have full control over which parameters are fitted via fit flags.

### Mathematical Framework

**Timing Model**:
```
t_model = t_pulsar + Σ(timing_params) + Σ(jump_params) + Σ(fdjump_params)
```

**Least-Squares Minimization**:
```
χ² = Σ[(t_obs - t_model)² / σ²]
```

**Design Matrix Construction**:
```
∂(residual)/∂(jump_param) = -1.0  (for JUMP)
∂(residual)/∂(fdjump_param) = -(f/f_ref)^index  (for FDJUMP)
∂(residual)/∂(fdjumpdm_param) = -(f/f_ref)^(-2) / DM_CONST  (for FDJUMPDM)
```

### Code Implementation

**Parameter Storage** (tempo2.h):
```c
// JUMP parameters
double jumpVal[MAX_JUMPS];        // Jump values
double jumpValErr[MAX_JUMPS];     // Jump uncertainties
int fitJump[MAX_JUMPS];           // Fit flags (1=fit, 0=fix)
char jumpStr[MAX_JUMPS][MAX_STRLEN]; // Jump descriptions

// FDJUMP parameters  
double fdjumpVal[MAX_JUMPS];      // FD jump values
double fdjumpValErr[MAX_JUMPS];   // FD jump uncertainties
int fdjumpIdx[MAX_JUMPS];         // Frequency indices
int fitfdJump[MAX_JUMPS];         // Fit flags
char fdjump_log;                  // Linear vs log scaling
```

**Fitting Functions** (t2fit_stdFitFuncs.C):
```c
// JUMP fitting function
double t2FitFunc_jump(pulsar *psr, int ipsr, double x, int ipos, param_label label, int k) {
    for (int l=0; l<psr[ipsr].obsn[ipos].obsNjump; l++) {
        if (psr[ipsr].obsn[ipos].jump[l]==k) {
            return -psr[ipsr].obsn[ipos].jumpScale[l];  // Phase jump
        }
    }
    return 0;
}

// FDJUMP fitting function
double t2FitFunc_fdjump(pulsar *psr, int ipsr, double x, int ipos, param_label label, int k) {
    for (int l=0; l<psr[ipsr].obsn[ipos].obsNfdjump; l++) {
        if (psr[ipsr].obsn[ipos].fdjump[l]==k) {
            int idx = psr[ipsr].fdjumpIdx[k];    
            if (idx == -2) {
                // DM jump: scales as ν^-2
                return pow(psr[ipsr].obsn[ipos].freqSSB/1e6, -2) / DM_CONST;
            } else {
                // Power-law frequency dependence
                return pow(psr[ipsr].obsn[ipos].freqSSB/1e9, idx);
            }
        }
    }
    return 0;
}
```

**Parameter Registration** (t2fit.C):
```c
// Register JUMP parameters for fitting
for (k=0; k<psr[p].nJumps; k++) {
    if (psr[p].fitJump[k] == 1) {
        registerFitFunc("JUMP", t2FitFunc_jump, psr[p].nJumps, k, "");
    }
}

// Register FDJUMP parameters for fitting  
for (k=0; k<psr[p].nfdJumps; k++) {
    if (psr[p].fitfdJump[k] == 1) {
        registerFitFunc("FDJUMP", t2FitFunc_fdjump, psr[p].nfdJumps, k, "");
    }
}
```

### User Control in Tempo2

**Full Parameter Control**:
- Users specify fit flags in .par file (1=fit, 0=fix)
- Parameters are updated iteratively during fitting
- Final fitted values and uncertainties written to .par file
- Standard least-squares uncertainty propagation

**Example Workflow**:
```bash
# 1. Set up .par file with JUMP parameters
echo "JUMP -sys system1 0.0 1" >> pulsar.par

# 2. Run Tempo2 fitting
tempo2 -f pulsar.par pulsar.tim

# 3. Updated .par file contains fitted JUMP values and uncertainties
```

## TempoNest Implementation

### Conceptual Approach
**Philosophy**: "Marginalize over nuisance parameters to focus on science parameters"

TempoNest treats JUMP-type parameters as nuisance parameters that should be analytically marginalized rather than sampled, allowing the analysis to focus on the parameters of scientific interest (red noise, DM variations, etc.).

### Mathematical Framework

**Bayesian Model**:
```
P(θ_science, θ_jumps | data) ∝ L(data | θ_science, θ_jumps) × P(θ_science) × P(θ_jumps)
```

**Analytical Marginalization**:
```
P(θ_science | data) = ∫ P(θ_science, θ_jumps | data) dθ_jumps
```

**JUMP Prior Distribution**:
```
θ_jump ~ N(μ_tempo2, σ_tempo2²)
```
Where μ and σ are taken from Tempo2 fitted values and uncertainties.

**Marginalized Likelihood**:
The JUMP parameters are integrated out analytically using their Gaussian priors, resulting in a reduced-dimensional likelihood that depends only on the science parameters.

### Code Implementation

**Parameter Detection** (timing_model.cpp):
```cpp
// Automatic detection of all fitted parameters
for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
    param_label p = fitinfo->paramIndex[iparam];
    const int k = fitinfo->paramCounters[iparam];
    
    if (p == param_JUMP) {
        t2_fitted_labels.push_back("JUMP_" + std::to_string(k));
        // Store fitted values and uncertainties for marginalization
        double jump_mean = globals::pulsar->jumpVal[k];
        double jump_err = globals::pulsar->jumpValErr[k] / error_scaling;
        t2_fit_values.push_back(jump_mean);
        t2_fit_errors.push_back(jump_err);
    } else {
        // Handle FDJUMP, FDJUMPDM, and other parameters generically
        t2_fitted_labels.push_back(globals::pulsar->param[p].shortlabel[k]);
        long double mean = globals::pulsar->param[p].prefit[k];
        long double err = globals::pulsar->param[p].err[k] / error_scaling;
        t2_fit_values.push_back(mean);
        t2_fit_errors.push_back(err);
    }
    
    // All parameters are marginalized by default
    marginalised.push_back(true);
}
```

**User Restriction Enforcement** (timing_model.cpp):
```cpp
void timing_model_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t&) {
    // Block all JUMP-type parameters from user control
    if (name.find("JUMP") != std::string::npos) {
        throw std::runtime_error("Cannot fit for jumps directly");
    }
    // ... rest of parameter setting logic
}
```

**Design Philosophy Display** (timing_model.h):
```cpp
void print_formatted(int element_index) const override {
    // Count JUMP-type parameters
    int jump_count = 0;
    for (size_t i = 0; i < marginalised.size(); i++) {
        if (marginalised[i] && t2_fitted_labels[i].find("JUMP") != std::string::npos) {
            jump_count++;
        }
    }
    
    // Display design philosophy
    if (jump_count > 0) {
        std::cout << "    └─ Design Philosophy: JUMP-type parameters (" << jump_count 
                  << ") are always marginalized, never sampled directly" << std::endl;
    }
}
```

### Marginalization Process

**Step 1**: TempoNest reads fitted JUMP-type parameters from Tempo2 output
**Step 2**: Parameters are assigned Gaussian priors based on fitted values/uncertainties  
**Step 3**: Analytical marginalization is performed in the likelihood calculation
**Step 4**: Users see only the science parameters in the sampling output

## Comparison: Tempo2 vs TempoNest

| Aspect | Tempo2 | TempoNest |
|--------|---------|-----------|
| **Treatment** | Normal fitted parameters | Nuisance parameters (marginalized) |
| **User Control** | Full control via fit flags | No user control (automatic) |
| **Mathematical Approach** | Least-squares fitting | Bayesian marginalization |
| **Parameter Updates** | Iterative optimization | Fixed at Tempo2 values |
| **Uncertainty Propagation** | Standard errors | Bayesian integration |
| **Output** | Fitted values + uncertainties | Marginalized out (not in chains) |
| **Design Philosophy** | "Fit what user specifies" | "Marginalize nuisances" |
| **Computational Cost** | O(N) per parameter | O(1) (analytical) |
| **Scientific Focus** | Complete parameter set | Science parameters only |

## Legacy TempoNest vs Modern TempoNest

| Feature | Legacy TempoNest | Modern TempoNest |
|---------|------------------|------------------|
| **JUMP Support** | ✅ Via doJumpMargin flag | ✅ Automatic marginalization |
| **FDJUMP Support** | ❌ Not implemented | ✅ Automatic marginalization |  
| **FDJUMPDM Support** | ❌ Not implemented | ✅ Automatic marginalization |
| **User Control** | doJumpMargin=0/1 toggle | No user control |
| **Parameter Detection** | Manual scanning | fitinfo structure |
| **Robustness** | Crash on complex .par files | Robust (after Option 3) |

## Practical Implications

### For Users

**When using Tempo2**:
- Set fit flags to control which JUMP parameters are fitted
- Monitor convergence and parameter correlations
- JUMP parameters appear in output .par file with uncertainties

**When using TempoNest**:
- JUMP-type parameters are automatically detected and marginalized
- No user configuration required or allowed
- Analysis focuses on science parameters (red noise, etc.)
- Faster sampling due to reduced parameter space

### For Developers

**Tempo2 Extension**:
- Add new JUMP-type parameters by extending the fitting framework
- Register new fit functions for parameter derivatives
- Update parameter storage structures

**TempoNest Extension**:
- New JUMP-type parameters are automatically supported
- No code changes needed if Tempo2 handles the parameter
- Focus development on science model elements

## Best Practices

### Choosing Between Tempo2 and TempoNest

**Use Tempo2 when**:
- You need precise JUMP parameter estimates
- Studying systematic offsets between observing systems
- Parameter correlations are scientifically important
- Building timing models for new pulsars

**Use TempoNest when**:
- Studying stochastic processes (red noise, DM variations)
- Gravitational wave analysis requiring noise characterization
- Bayesian model selection between noise models
- Uncertainty quantification on science parameters

### Recommended Workflow

1. **Initial Timing Model**: Use Tempo2 to establish timing parameters and JUMP values
2. **Noise Analysis**: Use TempoNest for Bayesian noise characterization
3. **Parameter Studies**: Return to Tempo2 for detailed parameter correlation studies
4. **Publication Results**: Combine both approaches for comprehensive analysis

## Technical Notes

### Error Scaling in TempoNest
```cpp
long double error_scaling = std::sqrt(globals::pulsar->fitChisq / globals::pulsar->fitNfree);
```
TempoNest applies post-fit error scaling to account for model inadequacy, ensuring realistic parameter uncertainties for marginalization.

### Parameter Index Safety
The current TempoNest implementation is vulnerable to parameter indexing errors when complex .par files contain many JUMP-type parameters. **Option 3** (modern C++ parameter scanner) addresses this robustness issue.

### Future Enhancements
- **Selective marginalization**: Allow users to choose which JUMP parameters to marginalize vs sample
- **Hierarchical priors**: Model systematic relationships between JUMP parameters
- **Correlated JUMP priors**: Account for correlations between different JUMP types

## References

- **Tempo2**: Hobbs, Edwards & Manchester (2006), MNRAS, 369, 655
- **TempoNest**: Lentati et al. (2014), MNRAS, 437, 3004  
- **Bayesian Pulsar Timing**: van Haasteren & Vallisneri (2014), Phys. Rev. D, 90, 104012
- **Pulsar Timing Arrays**: Burke-Spolaor et al. (2019), A&ARv, 27, 5

---

*This documentation reflects the design decisions and implementation details of TempoNest as of the current codebase. For questions or clarifications, refer to the TempoNest development team.*