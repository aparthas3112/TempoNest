# Enhanced TempoNest Corner Plot Tool

A powerful, user-friendly tool for creating corner plots from TempoNest MultiNest chains output. The enhanced version provides automatic parameter detection, professional formatting, and flexible parameter selection.

## Key Improvements

### 🚀 **Ease of Use**
- **One-command plotting**: `python plot_corner.py chains.txt` 
- **Automatic parameter detection**: No need to manually specify parameter names
- **Smart defaults**: Automatically excludes likelihood column for large parameter sets
- **Comprehensive help**: Built-in examples and documentation

### 🏷️ **Professional Formatting**
- **LaTeX labels**: Beautiful mathematical notation (e.g., $\log_{10} A_{red}$, $\gamma_{red}$)
- **Standard TempoNest naming**: Recognizes common parameter patterns
- **Automatic transformations**: EFAC converted from log₁₀ to linear scale
- **Publication-ready**: High-quality PDF/PNG output

### 🎯 **Flexible Parameter Selection**
- **By index**: `--params 0,1,2` 
- **By name**: `--params red_amp,red_gamma,efac`
- **Mixed selection**: Combine indices and names
- **Smart validation**: Clear error messages for invalid selections

### 📊 **Enhanced Statistics**
- **Quantiles**: Shows 16th, 50th, 84th percentiles
- **Smooth histograms**: Better visual appearance
- **Truth values**: Easy comparison with known values
- **Custom titles**: Add analysis context

## Installation Requirements

```bash
pip install numpy matplotlib corner
```

## Usage Examples

### Basic Usage
```bash
# Plot all parameters with auto-detection
python plot_corner.py test_mnest_chains.txt

# Output: Plots all 18 parameters (excludes likelihood automatically)
```

### Parameter Selection
```bash
# Plot first 3 parameters by index
python plot_corner.py chains.txt --params 0,1,2

# Plot red noise parameters by name  
python plot_corner.py chains.txt --params red_amp,red_gamma,efac

# Mix indices and names
python plot_corner.py chains.txt --params 0,red_gamma,17
```

### Truth Values
```bash
# Include known truth values for comparison
python plot_corner.py chains.txt --params 0,1,2 --truth -15.0,3.5,0.1

# Partial truth values (rest will be None)
python plot_corner.py chains.txt --params 0,1,2,3 --truth -15.0,3.5
```

### Output Options
```bash
# Save to file instead of displaying
python plot_corner.py chains.txt --output corner_plot.png

# High-resolution PDF for publication
python plot_corner.py chains.txt --output figure.pdf --dpi 300

# Add custom title
python plot_corner.py chains.txt --title "PSR J1234+5678 Red Noise Analysis"
```

### Advanced Options
```bash
# Force include likelihood column
python plot_corner.py chains.txt --include-likelihood

# Disable automatic transformations
python plot_corner.py chains.txt --no-transform

# Complete analysis with all options
python plot_corner.py chains.txt \
    --params red_amp,red_gamma,efac \
    --truth -15.0,3.5,0.1 \
    --title "Red Noise + EFAC Analysis" \
    --output analysis_corner.pdf \
    --dpi 300
```

## Automatic Parameter Detection

The tool automatically recognizes common TempoNest parameter patterns:

| Pattern | Parameters | Use Case |
|---------|------------|----------|
| **Red Noise Only** | `red_amp`, `red_gamma`, `efac` | Simple noise analysis |
| **Red + DM Noise** | `red_amp`, `red_gamma`, `dm_amp`, `dm_gamma`, `efac`, `equad` | Dual noise processes |
| **Full Model** | Noise + timing parameters + `efac` | Complete pulsar analysis |

### Recognized Parameter Names

**Noise Parameters:**
- `red_amp`, `red_amplitude`, `redamp` → $\log_{10} A_{red}$
- `red_gamma`, `red_spectral_index`, `redgamma` → $\gamma_{red}$
- `dm_amp`, `dm_amplitude`, `dmamp` → $\log_{10} A_{DM}$
- `dm_gamma`, `dm_spectral_index`, `dmgamma` → $\gamma_{DM}$

**Error Model:**
- `efac`, `efac_global` → EFAC (auto-converts from log to linear)
- `equad`, `equad_global` → EQUAD

**Timing Parameters:**
- `f0`, `f1`, `f2` → Frequency derivatives
- `raj`, `decj` → Position
- `pmra`, `pmdec` → Proper motion
- `pb`, `a1`, `e`, `omega` → Orbital elements

## File Format Support

The tool works with standard MultiNest output files:
- `chains.txt` - Raw chains
- `post_equal_weights.dat` - Equal-weighted posterior
- Any whitespace-delimited numerical file

**Expected format**: Each row is a sample, each column is a parameter, with likelihood typically in the last column.

## Error Handling

The enhanced tool provides clear error messages:
- **Missing files**: File not found warnings
- **Invalid parameters**: Unknown parameter names or indices
- **Format issues**: Data loading problems
- **Missing dependencies**: Package installation guidance

## Migration from Original Script

**Old workflow:**
```bash
python plot_corner.py \
    --input chains.txt \
    --params "red_amp,red_gamma,efac,likelihood" \
    --truth "-15.0,3.5,0.1"
```

**New workflow:**
```bash
python plot_corner.py chains.txt \
    --params red_amp,red_gamma,efac \
    --truth -15.0,3.5,0.1
```

**Key differences:**
- ✅ Input file is positional argument (simpler)
- ✅ No quotes needed around parameter lists  
- ✅ Automatic parameter detection
- ✅ Likelihood excluded by default
- ✅ Better error handling and validation

## Tips for Best Results

1. **Parameter Selection**: Use descriptive names when possible for clearer plots
2. **Truth Values**: Always include if available for validation
3. **Output Format**: Use PDF for publications, PNG for presentations
4. **Large Models**: Let auto-detection exclude likelihood for cleaner plots
5. **Custom Labels**: The tool handles LaTeX formatting automatically

## Troubleshooting

**"Parameter not found" warnings:**
- Check parameter spelling and case sensitivity
- Use `--params` without arguments to see auto-detected names
- Fall back to indices if names don't match

**"corner package not found":**
```bash
pip install corner
```

**Empty or malformed plots:**
- Verify file format (whitespace-delimited numbers)
- Check for NaN or infinite values in chains
- Ensure sufficient samples for meaningful statistics

The enhanced corner plot tool transforms TempoNest analysis visualization from a cumbersome process into a simple, one-command operation while providing professional-quality output suitable for publications.