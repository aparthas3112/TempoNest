# TempoNest Python Tools

## Overview

TempoNest includes Python-based tools for configuration generation and result analysis. These tools provide user-friendly interfaces for creating JSON configurations and visualizing analysis results.

## Table of Contents

1. [Configuration GUI (TempoNest_JSON.py)](#configuration-gui)
2. [Corner Plot Tool (plot_corner.py)](#corner-plot-tool)
3. [Installation and Setup](#installation-and-setup)
4. [Usage Examples](#usage-examples)

---

## Configuration GUI

The TempoNest Configuration GUI is a Streamlit-based web interface for creating and validating JSON configuration files.

### Features

**Core Functionality:**
- Interactive JSON configuration generation
- Real-time validation and error checking
- Export to file or copy to clipboard
- Parameter range suggestions and validation

**Model Support:**
- **Three timing model modes**: Marginalize (default), Fit All, Manual
- **All model elements**: Red noise, DM noise, EFAC, EQUAD, ECORR, Solar wind
- **Enhanced ECORR support**: Multi-flag selection with shared priors
- **Element type validation**: Parameter ranges reset correctly when changing types

**User Experience:**
- **Improved defaults**: Per-flag EFAC/EQUAD, 500 live points, local output directory
- **Parameter tooltips**: Helpful descriptions and suggested ranges
- **Enhanced solar wind support**: Proper default ranges and documentation
- **Configuration templates**: Quick start with common analysis types

### Starting the GUI

```bash
# From the TempoNest directory
cd src/python/temponest
streamlit run TempoNest_JSON.py

# Alternative: Direct path
streamlit run src/python/temponest/TempoNest_JSON.py
```

The GUI will open in your default web browser at `http://localhost:8501`.

### GUI Sections

#### 1. Global Settings
- **Use Original Errors**: Toggle between original and Tempo2-fitted uncertainties
- **Tempo2 Iterations**: Number of pre-analysis fitting iterations
- **Test Mode**: Enable debugging output

#### 2. Sampler Configuration
- **Live Points**: Choose from preset values (500/1000/4000) or custom
- **Efficiency**: Exploration vs speed trade-off
- **Output Settings**: File prefix and directory
- **Advanced Options**: Evidence tolerance, resume capability

#### 3. Model Elements
- **Add Elements**: Dropdown selection of all supported model types
- **Parameter Configuration**: Interactive parameter setup with validation
- **Element Ordering**: Drag-and-drop reordering of elements
- **Validation**: Real-time checking of parameter ranges and compatibility

#### 4. Export Options
- **Download JSON**: Save configuration to file
- **Copy to Clipboard**: Quick copying for command-line use
- **Validation Report**: Check configuration before export

### Advanced Features

#### ECORR Multi-Flag Configuration

The GUI supports advanced ECORR setup:

1. **Flag Selection**: Choose from available flag types (-B, -fe, -f, -sys, -chan, -group)
2. **Shared Priors**: Apply same prior ranges across selected flags
3. **Automatic Backend Detection**: System discovers backends within each flag
4. **Legacy Compatibility**: Default settings match legacy TempoNest behavior

#### Parameter Range Validation

- **Amplitude Ranges**: Suggested ranges based on typical pulsar values
- **Spectral Index Limits**: Physical bounds for red and DM noise
- **Error Scaling**: Reasonable EFAC/EQUAD ranges for instrumentation
- **Cross-Validation**: Check for conflicting or redundant elements

---

## Corner Plot Tool

Enhanced corner plot generation with multi-dataset comparison capabilities.

### Basic Usage

```bash
# Simple corner plot with automatic parameter detection
python plot_corner.py chains.txt

# Specify output file
python plot_corner.py chains.txt --output corner_plot.png

# Custom title and parameter subset
python plot_corner.py chains.txt --title "Red Noise Analysis" --params 0,1,2
```

### Advanced Comparison Mode

Compare posteriors from multiple analyses:

#### Two-Dataset Comparison

```bash
# Compare two analyses with same parameter structure
python plot_corner.py chains1.txt --compare /path/to/dir2/ \
    --param-map-1 red_amp,red_gamma,efac \
    --param-map-2 red_amp,red_gamma,efac
```

#### Different Parameter Orderings

```bash
# Handle different parameter orders between analyses
python plot_corner.py chains1.txt --compare /path/to/dir2/ \
    --param-map-1 0,1,2 \
    --param-map-2 2,0,1
```

#### Multiple Comparisons

```bash
# Compare up to 10 datasets simultaneously
python plot_corner.py main_chains.txt \
    --compare dir2/ --compare dir3/ \
    --param-map-1 red_amp,red_gamma,efac \
    --param-map-2 red_amp,red_gamma,efac \
    --param-map-2 0,1,2 \
    --label-1 "GPU Analysis" \
    --label-2 "CPU Test" \
    --label-2 "Legacy"
```

### Features

**Parameter Recognition:**
- **Automatic Detection**: Recognizes common parameter names and applies LaTeX formatting
- **Solar Wind Support**: Enhanced recognition for solar wind parameters
- **Custom Labels**: User-defined parameter names and units
- **Missing Parameter Handling**: Graceful handling of incomplete parameter files

**Visual Features:**
- **10-Color Palette**: Distinct colors for multi-dataset comparison
- **Transparency Control**: Adjustable alpha for overlapping distributions
- **Flexible Layouts**: Automatic subplot arrangement for different parameter counts
- **High-Quality Output**: Publication-ready figures with customizable DPI

**File Handling:**
- **Multiple Formats**: Support for various chain file formats
- **Parameter File Detection**: Automatically finds corresponding parameter name files
- **Directory Structure**: Intelligent file discovery in analysis directories
- **Error Handling**: Robust handling of missing or malformed files

### Command-Line Options

```bash
python plot_corner.py [PRIMARY_CHAINS] [OPTIONS]

Required:
  PRIMARY_CHAINS        Main chains file to plot

Optional:
  --compare PATH        Additional dataset directory (can be repeated)
  --param-map-1 LIST    Parameter mapping for primary dataset
  --param-map-2 LIST    Parameter mapping for comparison dataset (can be repeated)
  --label-1 TEXT        Label for primary dataset
  --label-2 TEXT        Label for comparison dataset (can be repeated)
  --output PATH         Output file path (default: auto-generated)
  --title TEXT          Plot title
  --params INDICES      Subset of parameters to plot (e.g., "0,1,2")
  --dpi INTEGER         Output resolution (default: 300)
  --alpha FLOAT         Transparency for overlapping distributions (default: 0.7)
```

### Parameter Mapping

Parameter mapping supports both names and indices:

**By Name:**
```bash
--param-map-1 "red_amplitude,red_spectral_index,efac_group"
```

**By Index:**
```bash
--param-map-1 "0,1,2"
```

**Mixed (discouraged but supported):**
```bash
--param-map-1 "red_amplitude,1,efac_group"
```

### File Structure Requirements

The tool expects this directory structure:

```
analysis_directory/
├── chains.txt                    # Posterior samples (required)
├── post_equal_weights.dat        # Alternative chain file name
├── paramnames.txt                # Parameter names (optional)
├── TNest-paramnames.txt          # TempoNest parameter names (optional)
└── stats.dat                     # Analysis statistics (optional)
```

---

## Installation and Setup

### Dependencies

Install required Python packages:

```bash
# Essential packages
pip install streamlit pandas numpy matplotlib corner

# Optional: Enhanced plotting
pip install seaborn plotly

# For development
pip install jupyter notebook ipython
```

### Streamlit Configuration

Create `~/.streamlit/config.toml` for optimal experience:

```toml
[server]
port = 8501
headless = false

[browser]
gatherUsageStats = false

[theme]
primaryColor = "#FF6B6B"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
```

### Environment Setup

Add to your shell configuration:

```bash
# TempoNest Python tools
export PYTHONPATH="${PYTHONPATH}:/path/to/TempoNest/src/python"

# Streamlit configuration
export STREAMLIT_SERVER_HEADLESS=false
export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
```

---

## Usage Examples

### Complete Workflow Example

1. **Generate Configuration**:
   ```bash
   streamlit run src/python/temponest/TempoNest_JSON.py
   # Use GUI to create config.json
   ```

2. **Run Analysis**:
   ```bash
   tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json
   ```

3. **Create Corner Plot**:
   ```bash
   python src/python/temponest/plot_corner.py \
       TNest-post_equal_weights.dat \
       --title "Pulsar Noise Analysis" \
       --output analysis_corner.png
   ```

### Comparison Analysis Workflow

1. **Run Multiple Analyses**:
   ```bash
   # Analysis 1: GPU
   tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config_gpu.json
   mv TNest-* analysis_gpu/
   
   # Analysis 2: CPU
   tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config_cpu.json
   mv TNest-* analysis_cpu/
   ```

2. **Compare Results**:
   ```bash
   python src/python/temponest/plot_corner.py \
       analysis_gpu/TNest-post_equal_weights.dat \
       --compare analysis_cpu/ \
       --label-1 "GPU Analysis" \
       --label-2 "CPU Analysis" \
       --output gpu_vs_cpu_comparison.png
   ```

### Batch Processing Example

Create multiple configurations and analyses:

```bash
#!/bin/bash
# batch_analysis.sh

# Generate configurations for different noise models
python generate_configs.py  # Custom script using GUI components

# Run analyses
for config in configs/*.json; do
    output_dir="results/$(basename $config .json)"
    mkdir -p $output_dir
    
    tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile $config
    mv TNest-* $output_dir/
    
    # Generate corner plot
    python src/python/temponest/plot_corner.py \
        $output_dir/TNest-post_equal_weights.dat \
        --output $output_dir/corner_plot.png
done
```

## Troubleshooting

### Common Issues

**Streamlit won't start:**
```bash
# Check if port is in use
lsof -i :8501

# Use different port
streamlit run TempoNest_JSON.py --server.port 8502
```

**Import errors:**
```bash
# Check Python path
python -c "import sys; print(sys.path)"

# Install missing packages
pip install --user streamlit corner matplotlib
```

**Plot generation fails:**
```bash
# Check file permissions
ls -la *chains*

# Verify file format
head -5 TNest-post_equal_weights.dat
```

**Parameter mapping errors:**
```bash
# Check parameter file
cat TNest-paramnames.txt

# Use indices instead of names
--param-map-1 "0,1,2,3"
```

### Performance Tips

**Large Chain Files:**
- Use parameter subset (`--params`) for faster plotting
- Reduce DPI for preview plots
- Consider data thinning for very large chains

**Memory Usage:**
- Close unused browser tabs when using Streamlit
- Use `--alpha 0.5` for better performance with many datasets
- Monitor memory usage with large comparison plots

**Streamlit Performance:**
- Restart Streamlit app if it becomes slow
- Clear browser cache periodically
- Use local file uploads instead of network paths

See [troubleshooting.md](troubleshooting.md) for additional support and [examples.md](examples.md) for complete configuration examples.