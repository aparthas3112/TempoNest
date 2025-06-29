# TempoNest Troubleshooting Guide

## Overview

This guide covers common issues, error messages, and debugging procedures for TempoNest installation, configuration, and runtime problems.

## Table of Contents

1. [Build and Installation Issues](#build-and-installation-issues)
2. [Runtime Errors](#runtime-errors)
3. [Performance Issues](#performance-issues)
4. [GPU-Specific Problems](#gpu-specific-problems)
5. [Configuration and JSON Issues](#configuration-and-json-issues)
6. [Scientific Validation](#scientific-validation)
7. [Debug Mode and Logging](#debug-mode-and-logging)

---

## Build and Installation Issues

### Tempo2 Not Found

**Error Messages:**
```bash
configure: error: TEMPO2 not found
```

**Solutions:**
```bash
# Set environment variable
export TEMPO2=/path/to/tempo2
export PATH=$TEMPO2/bin:$PATH

# Verify installation
echo $TEMPO2
which tempo2
$TEMPO2/bin/tempo2 --help
```

**Common Locations:**
- `/usr/local/tempo2`
- `/opt/tempo2`
- `$HOME/software/tempo2`
- Module systems: `module load tempo2; echo $TEMPO2`

### Missing Dependencies

**GSL Not Found:**
```bash
# Ubuntu/Debian
sudo apt-get install libgsl-dev

# CentOS/RHEL
sudo yum install gsl-devel

# macOS
brew install gsl
```

**BLAS/LAPACK Issues:**
```bash
# Ubuntu/Debian
sudo apt-get install libopenblas-dev liblapack-dev

# CentOS/RHEL
sudo yum install openblas-devel lapack-devel

# Specify in configure
./configure --with-blas=openblas
```

**Eigen3 Headers Missing:**
```bash
# Ubuntu/Debian
sudo apt-get install libeigen3-dev

# Manual specification
./configure --with-eigen=/usr/include/eigen3
```

### MPI Compiler Issues

**Error Messages:**
```bash
configure: error: MPI compiler not found
```

**Solutions:**
```bash
# Install MPI
sudo apt-get install libopenmpi-dev openmpi-bin

# Set compiler explicitly
export CC=mpicc
export CXX=mpicxx

# Verify MPI installation
which mpicc
mpicc --version
```

### ArrayFire Configuration

**ArrayFire Not Found:**
```bash
# Option 1: Disable GPU support
./configure --without-arrayfire

# Option 2: Specify ArrayFire path
export ARRAYFIRE_PATH=/path/to/arrayfire
./configure --with-arrayfire=$ARRAYFIRE_PATH

# Option 3: System installation
sudo apt-get install arrayfire-dev  # if available
```

**ArrayFire Library Issues:**
```bash
# Check installation
ls $ARRAYFIRE_PATH/lib64/libaf*
export LD_LIBRARY_PATH=$ARRAYFIRE_PATH/lib64:$LD_LIBRARY_PATH

# Test ArrayFire
export AF_PRINT_ERRORS=1
```

### Autotools Problems

**Missing autogen.sh:**
```bash
# Regenerate autotools files
autoreconf -fiv

# Or manually
autoconf
automake --add-missing
```

**Libtool Issues:**
```bash
# Install libtool
sudo apt-get install libtool

# Check installation
which libtoolize
libtoolize --version
```

---

## Runtime Errors

### Plugin Loading Failures

**Error Messages:**
```bash
ERROR: Cannot load plugin temponest
```

**Solutions:**
```bash
# Reinstall plugin
make temponest-install

# Check plugin directory
ls $TEMPO2/plugins/temponest_*

# Verify plugin installation
echo $TEMPO2/plugins
$TEMPO2/bin/tempo2 -gr temponest --help
```

**Library Path Issues:**
```bash
# Add library paths
export LD_LIBRARY_PATH=$TEMPO2/lib:$LD_LIBRARY_PATH

# For ArrayFire
export LD_LIBRARY_PATH=$ARRAYFIRE_PATH/lib64:$LD_LIBRARY_PATH

# For MultiNest
export LD_LIBRARY_PATH=/path/to/multinest:$LD_LIBRARY_PATH
```

### Segmentation Faults

**Debug with GDB:**
```bash
# Compile with debug symbols
make clean
./configure --enable-debug
make temponest

# Run with debugger
gdb --args $TEMPO2/bin/tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json
(gdb) run
(gdb) bt  # backtrace when crash occurs
```

**Common Causes:**
1. **Malformed .par files**: Check for unusual parameter entries
2. **Invalid JSON**: Validate configuration file syntax
3. **Memory issues**: Monitor with `valgrind` or `AddressSanitizer`
4. **Library version conflicts**: Ensure consistent dependency versions

### File Access Issues

**Cannot Read Input Files:**
```bash
# Check file permissions
ls -la pulsar.par pulsar.tim config.json

# Verify file paths
realpath pulsar.par
head -5 pulsar.par
```

**Output Directory Problems:**
```bash
# Create output directory
mkdir -p results/

# Check write permissions
touch results/test_file && rm results/test_file

# Update config.json
"output_root": "results/TNest-"
```

---

## Performance Issues

### Slow Convergence

**MultiNest Parameters:**
```json
{
  "sampler": {
    "live_points": 1000,     // Increase for better exploration
    "efficiency": 0.1,       // Decrease for more thorough sampling
    "evidence_tolerance": 0.5 // Decrease for higher precision
  }
}
```

**Diagnostic Checks:**
```bash
# Monitor acceptance rate (should be 5-15%)
tail -f TNest-stats.dat

# Check parameter ranges
grep "min_value\|max_value" config.json
```

### High Memory Usage

**Memory Monitoring:**
```bash
# Monitor during run
top -p $(pgrep tempo2)
watch -n 5 'free -h'

# Check for memory leaks
valgrind --tool=memcheck --leak-check=full tempo2 -gr temponest ...
```

**Reduction Strategies:**
```json
{
  "elements": [
    {
      "name": "Power Law Red Noise",
      "days_per_coeff": 50.0  // Increase to reduce coefficients
    }
  ]
}
```

### CPU Utilization

**MPI Scaling:**
```bash
# Use multiple cores
mpirun -np 8 tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json

# Check MPI performance
export OMPI_DEBUG=1
```

**Thread Optimization:**
```bash
# OpenMP settings
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
```

---

## GPU-Specific Problems

### CUDA Driver Issues

**Check CUDA Installation:**
```bash
nvidia-smi
nvcc --version
cat /proc/driver/nvidia/version
```

**Driver Updates:**
```bash
# Ubuntu
sudo apt-get update
sudo apt-get install nvidia-driver-XXX

# Check compatibility
nvidia-smi | grep "Driver Version"
```

### ArrayFire GPU Problems

**Enable Debugging:**
```bash
export AF_DEBUG=1
export AF_PRINT_ERRORS=1
export AF_TIMER=1
```

**Backend Selection:**
```bash
# Force CUDA backend
export AF_DEFAULT_DEVICE=CUDA

# Check available backends
./arrayfire_test_backends
```

**Memory Issues:**
```bash
# Monitor GPU memory
watch -n 1 nvidia-smi

# Reduce problem size if GPU memory is insufficient
"days_per_coeff": 60.0  // Fewer coefficients
```

### GPU Performance Issues

**Low GPU Utilization:**
```bash
# Monitor GPU usage
nvidia-smi dmon -s pucvmet -i 0

# Profile GPU usage
nsys profile --stats=true tempo2 -gr temponest ...
```

**Common Causes:**
1. **Small problem size**: GPU overhead dominates for <20 parameters
2. **Memory bandwidth**: Check for CPU-GPU transfer bottlenecks
3. **ArrayFire configuration**: Verify optimal backend selection

---

## Configuration and JSON Issues

### JSON Syntax Errors

**Validation Tools:**
```bash
# Command-line validation
python -m json.tool config.json

# Online validators
# https://jsonlint.com/
```

**Common Syntax Issues:**
```json
// Incorrect (trailing comma)
{
  "live_points": 1000,
  "efficiency": 0.1,
}

// Correct
{
  "live_points": 1000,
  "efficiency": 0.1
}
```

### Parameter Range Issues

**Invalid Ranges:**
```json
// Error: min_value > max_value
{
  "name": "amplitude",
  "min_value": -10,
  "max_value": -18  // Should be > min_value
}
```

**Scientific Validation:**
```json
// Reasonable red noise amplitude ranges
{
  "name": "amplitude",
  "min_value": -18,    // Very weak noise
  "max_value": -10     // Very strong noise
}
```

### Element Configuration Errors

**Missing Required Fields:**
```json
// Error: Missing prior_type
{
  "name": "amplitude",
  "min_value": -18,
  "max_value": -10
  // "prior_type": "log_uniform"  // Required!
}
```

**Flag Validation:**
```bash
# Check available flags in data
grep "^[[:space:]]*[0-9]" pulsar.tim | cut -d' ' -f4- | sort | uniq

# Common flags: -group, -B, -fe, -f, -sys
```

---

## Scientific Validation

### Result Verification

**Cross-Check with Known Results:**
```bash
# Use reference test case
tempo2 -gr temponest -f tests/test_data/test.par tests/test_data/test.tim -Cfile tests/test_data/test.json

# Expected likelihood: 929.961754 ± 1e-6
grep "Log Likelihood" TNest-stats.dat
```

**Parameter Sanity Checks:**
```bash
# Check posterior distributions
python plot_corner.py TNest-post_equal_weights.dat

# Verify parameter ranges are reasonable
python -c "
import numpy as np
chains = np.loadtxt('TNest-post_equal_weights.dat')
print('Parameter ranges:')
for i in range(chains.shape[1]):
    print(f'  Param {i}: [{chains[:,i].min():.3f}, {chains[:,i].max():.3f}]')
"
```

### Reproducibility Testing

**Fixed Seed Runs:**
```json
{
  "sampler": {
    "seed": 12345,  // Fixed seed for reproducibility
    "live_points": 500
  }
}
```

**Statistical Consistency:**
```bash
# Run multiple times with different seeds
for seed in 1 2 3 4 5; do
    sed "s/\"seed\": [0-9]*/\"seed\": $seed/" config.json > config_$seed.json
    tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config_$seed.json
    mv TNest-* results_seed_$seed/
done

# Compare evidence values (should be within 2-3σ)
```

---

## Debug Mode and Logging

### Enable Debug Output

**Global Debug Settings:**
```json
{
  "globals": {
    "test_mode": true,     // Enable debug output
    "use_original_errors": true,
    "num_tempo2_its": 1
  }
}
```

**Verbose Sampler Output:**
```json
{
  "sampler": {
    "verbose": true,       // Enable MultiNest verbose output
    "sample": true
  }
}
```

### Log File Analysis

**MultiNest Log Files:**
```bash
# Monitor sampling progress
tail -f TNest-stats.dat

# Check final evidence
grep "Global Evidence" TNest-stats.dat

# Monitor live points
wc -l TNest-live.points
```

**Error Log Analysis:**
```bash
# Search for error patterns
grep -i "error\|warning\|fail" tempo2_output.log

# Check memory usage patterns
grep -i "memory\|malloc\|alloc" tempo2_output.log
```

### Custom Debugging

**Add Debug Prints:**
```cpp
// In source code for detailed debugging
logdbg("Parameter value: %f", parameter_value);
logdbg("Matrix dimensions: %d x %d", rows, cols);
```

**Conditional Compilation:**
```cpp
#ifdef DEBUG
    std::cout << "Debug: entering function" << std::endl;
    std::cout << "Debug: variable value = " << value << std::endl;
#endif
```

### Performance Profiling

**Timing Analysis:**
```bash
# Time the complete run
time tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json

# Profile with gprof (if compiled with -pg)
gprof tempo2 gmon.out > profile_report.txt
```

**Memory Profiling:**
```bash
# Valgrind memory analysis
valgrind --tool=massif tempo2 -gr temponest -f pulsar.par pulsar.tim -Cfile config.json

# Analyze memory usage
ms_print massif.out.XXXXX
```

## Getting Help

### Community Support

**Documentation:**
- Read all documentation files in `docs/` directory
- Check [examples.md](examples.md) for working configurations
- Review [models.md](models.md) for model-specific issues

**Issue Reporting:**
When reporting issues, include:
1. **System Information**: OS, compiler versions, dependencies
2. **Build Configuration**: Configure options and environment variables
3. **Error Messages**: Complete error output with context
4. **Reproducible Example**: Minimal case that demonstrates the problem
5. **Expected vs Actual**: What you expected vs what happened

**Debugging Checklist:**
- [ ] Environment variables set correctly (`$TEMPO2`, paths)
- [ ] All dependencies installed and accessible
- [ ] JSON configuration validated
- [ ] Input files readable and well-formed
- [ ] Output directory writable
- [ ] Debug mode enabled for detailed output
- [ ] Known test case works correctly

For model-specific questions, see [models.md](models.md). For GPU issues, see [gpu_optimization.md](gpu_optimization.md). For Python tools, see [python_tools.md](python_tools.md).