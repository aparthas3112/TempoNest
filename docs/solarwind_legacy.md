# Solar Wind Model Legacy Documentation - TempoNest

## Overview

TempoNest implements a two-component solar wind model based on the You et al. (2007) formalism:
1. **Deterministic Component**: Corrects for systematic solar wind variations
2. **Stochastic Component**: Models unaccounted solar wind fluctuations as additional white noise

Both components scale with the expected dispersion measure dependence (∝ 1/f²) and geometric factors.

## Mathematical Formulation

### Deterministic Solar Wind Correction
Applied directly to timing residuals:

```
δt_det(t) = (Δn_e - n_e,ref) × τ_dis²(t)
```

Where:
- `Δn_e` = fitted solar wind electron density parameter (`SolarWind`)
- `n_e,ref` = reference solar wind electron density (`pulse->ne_sw`)  
- `τ_dis²(t)` = dispersion delay squared term (`pulse->obsn[o].tdis2`)

### Stochastic Solar Wind Noise
Added to white noise covariance:

```
σ²_sw(t) = [A_sw × τ_dis²(t) / n_e,ref]²
```

Where:
- `A_sw` = solar wind noise amplitude (`WhiteSolarWind = 10^(fitted_param)`)

## Software Architecture

### Configuration System

**Setup Parameters** (from `TempoNestPlugin.cpp`):
```cpp
int FitSolarWind = 0;           // Enable deterministic component
int FitWhiteSolarWind = 0;      // Enable stochastic component  
double *SolarWindPrior;         // Prior bounds for deterministic param
double *WhiteSolarWindPrior;    // Prior bounds for stochastic param
```

**Configuration Integration**:
- Read from config file via `setupparams()`
- Stored in `MNStruct` context
- Controls parameter counting for sampling dimensions

### Data Structures

**Key Variables in MNStruct**:
```cpp
struct MNStruct {
    int FitSolarWind;           // Flag: fit deterministic component
    int FitWhiteSolarWind;      // Flag: fit stochastic component
    pulsar *pulse;              // Contains ne_sw, tdis2 per observation
    // ... other fields
};
```

**Per-Observation Data** (in `pulse->obsn[o]`):
```cpp
double ne_sw;                   // Reference solar wind density
double tdis2;                   // Dispersion delay squared term
```

## Implementation Details

### Parameter Extraction
Located in `NewLRedMarginLogLike()` around line 680:

```cpp
double SolarWind=0;
double WhiteSolarWind = 0;

if(((MNStruct *)globalcontext)->FitSolarWind == 1){
    SolarWind = Cube[pcount];           // Linear parameter
    pcount++;
}
if(((MNStruct *)globalcontext)->FitWhiteSolarWind == 1){
    WhiteSolarWind = pow(10.0, Cube[pcount]);  // Log10 parameter
    pcount++;
}
```

### Noise Model Integration
In white noise calculation (around line 770):

```cpp
double SWTerm = WhiteSolarWind*((MNStruct *)globalcontext)->pulse->obsn[o].tdis2/
                ((MNStruct *)globalcontext)->pulse->ne_sw;

// Added to inverse noise variance
Noise[o]= 1.0/(pow(EFACterm,2) + EQUAD[...] + ShannonJitterTerm + 
               SWTerm*SWTerm + DMEQUADTerm*DMEQUADTerm);
```

### Deterministic Correction
Applied to residuals (around line 1250):

```cpp
if(((MNStruct *)globalcontext)->FitSolarWind == 1){
    for(int o=0;o<((MNStruct *)globalcontext)->pulse->nobs; o++){
        Resvec[o]-= (SolarWind-((MNStruct *)globalcontext)->pulse->ne_sw)*
                    ((MNStruct *)globalcontext)->pulse->obsn[o].tdis2;
    }
}
```

## Data Flow Architecture

### Initialization Phase
1. **Config Loading**: `setupparams()` reads solar wind flags from config file
2. **Parameter Setup**: Prior bounds set in `Dpriors` array  
3. **Context Assignment**: Flags stored in `MNStruct` for likelihood access
4. **Dimension Counting**: `swdims` incremented for each enabled component

### Likelihood Evaluation Phase
1. **Parameter Extraction**: Extract `SolarWind`/`WhiteSolarWind` from parameter cube
2. **Deterministic Correction**: Subtract systematic delays from residuals
3. **Noise Calculation**: Add stochastic component to noise covariance per observation
4. **Likelihood Computation**: Use modified residuals and noise in Gaussian likelihood

### Parameter Flow
```
Config File → setupparams() → MNStruct → NewLRedMarginLogLike()
     ↓
Prior Setup → Dpriors → Sampler → Parameter Cube → Solar Wind Components
```

## Integration Points

### With Tempo2
- **Data Source**: `pulse->obsn[o].tdis2`, `pulse->ne_sw` from Tempo2 infrastructure
- **Units**: Assumes Tempo2's dispersion delay conventions
- **Coordinates**: Uses Tempo2's geometric calculations for solar wind path

### With Noise Model
- **White Noise Stack**: Stochastic component added to total white noise variance
- **Independence**: Assumes independence from other noise components  
- **Scaling**: Multiplicative with existing error budget

### With Sampling
- **Parameter Space**: Linear for deterministic, log10 for stochastic
- **Prior Integration**: Uses standard TempoNest prior system
- **Marginalization**: No special marginalization handling

## Critical Implementation Notes

### Numerical Considerations
- **Log Parameters**: `WhiteSolarWind` uses log10 parameterization for numerical stability
- **Zero Handling**: No explicit zero-protection for `ne_sw` division
- **Units**: Implicit assumption about `tdis2` units matching delay expectations

### Physical Assumptions
- **Linear Scaling**: Deterministic component scales linearly with density difference
- **Frequency Dependence**: Encoded in `tdis2` calculation (should be ∝ 1/f²)
- **Geometric Scaling**: Solar wind path effects embedded in `tdis2`

### Memory Layout
- **Per-Observation**: Both `tdis2` and `ne_sw` stored per observation
- **Global Parameters**: `SolarWind`/`WhiteSolarWind` are global fit parameters
- **No Caching**: Recalculated each likelihood evaluation

## GPU Implementation Considerations

### Parallelization Opportunities
- **Per-Observation**: Both deterministic correction and noise calculation are embarrassingly parallel
- **Memory Access**: Requires `tdis2`, `ne_sw` arrays accessible to GPU kernels
- **Parameter Broadcasting**: Global solar wind parameters need broadcast to all threads

### Memory Patterns
- **Observation Arrays**: `tdis2[nobs]`, `ne_sw[nobs]` need GPU memory allocation
- **Residual Modification**: In-place modification of `Resvec[nobs]` array
- **Noise Computation**: Component of per-observation noise calculation

### Computational Complexity
- **O(nobs)**: Linear scaling with number of observations
- **Minimal Computation**: Simple arithmetic operations per observation
- **No Dependencies**: Each observation computed independently

## Extension Points for Improved Models

### Current Limitations
- **Static Reference**: Single `ne_sw` reference value per observation
- **No Time Dependence**: No explicit modeling of solar wind time evolution  
- **Simple Scaling**: Linear relationship may be oversimplified
- **No Frequency Evolution**: `tdis2` pre-computed, no dynamic frequency dependence

### Potential Enhancements
- **Time-Variable Models**: Solar cycle, seasonal variations
- **Non-Linear Scaling**: More sophisticated density-delay relationships
- **Correlated Noise**: Temporal correlations in solar wind fluctuations
- **Multi-Component**: Separate fast/slow solar wind components

### Architecture Flexibility
- **Modular Design**: Current implementation easily replaceable
- **Parameter Extension**: Additional parameters can be added to parameter cube
- **Data Integration**: Additional per-observation data can be incorporated
- **GPU Kernels**: New models can leverage same parallelization patterns