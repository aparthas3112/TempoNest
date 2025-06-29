# TempoNest Installation Guide

## Overview

TempoNest is a Bayesian pulsar timing analysis tool that integrates with Tempo2 as a plugin. This guide covers installation, dependencies, and initial testing.

## Dependencies

### Required Dependencies

Before building TempoNest, ensure these components are installed:

**Core Requirements:**
- **Tempo2**: Set `$TEMPO2` environment variable to installation path
- **MultiNest**: Nested sampling library with Fortran bindings
- **GSL**: GNU Scientific Library
- **BLAS/LAPACK**: Linear algebra libraries (OpenBLAS recommended)
- **Eigen3**: C++ template library for linear algebra
- **MPI**: Message Passing Interface (mpicc/mpicxx compilers)

**Optional Dependencies:**
- **ArrayFire**: GPU acceleration library (for CUDA/OpenCL support)

### Environment Variables

Set the following environment variables before building:

```bash
# Required: Tempo2 installation path
export TEMPO2=/path/to/tempo2

# Optional: ArrayFire for GPU acceleration
export ARRAYFIRE_PATH=/path/to/arrayfire

# Optional: Custom library paths
export MULTINEST_PATH=/path/to/multinest
```

## Build Process

### 1. Initial Setup

Run this once after cloning the repository:

```bash
# Generate build system
./autogen.sh

# Configure build (basic)
./configure
```

### 2. Configure Options

Key configuration flags for customization:

```bash
# Basic configuration
./configure --with-blas=openblas

# Enable GPU support
./configure --with-blas=openblas --with-arrayfire=/path/to/arrayfire

# Debug build
./configure --enable-debug

# Custom installation paths
./configure --with-tempo2-plug-dir=/custom/path \
           --with-eigen=/path/to/eigen \
           --with-multinest=/path/to/multinest
```

**Common Configure Options:**
- `--enable-debug` - Debug build with symbols and verbose output
- `--with-eigen=/path` - Specify Eigen headers location
- `--with-arrayfire=/path` - Enable GPU support with ArrayFire
- `--with-tempo2-plug-dir=/path` - Custom plugin installation directory
- `--with-blas=openblas` - Use OpenBLAS for linear algebra

### 3. Compilation

```bash
# Build the plugin
make temponest

# Install plugin to Tempo2
make temponest-install

# Clean build artifacts
make clean
```

### 4. Automated Build Script

For convenience, use the provided build script:

```bash
# Basic build
./scripts/build_temponest.sh

# Clean build
./scripts/build_temponest.sh --clean

# Skip tests and installation
./scripts/build_temponest.sh --no-tests --no-install
```

## Testing

### Basic Testing

```bash
# Run all tests (builds first if needed)
make test

# Manual test with Tempo2
$TEMPO2/bin/tempo2 -gr temponest \
    -f tests/test_data/test.par \
    tests/test_data/test.tim \
    -cfile tests/test_data/test.json
```

### Verify Installation

Check that the plugin was installed correctly:

```bash
# List installed plugins
ls $TEMPO2/plugins/temponest_*_plug.t2

# Test plugin loading
echo "quit" | timeout 10s $TEMPO2/bin/tempo2 -gr temponest
```

## System Requirements

### Minimum Requirements

- **CPU**: Multi-core processor (4+ cores recommended)
- **RAM**: 8GB minimum (16GB+ for large datasets)
- **Storage**: 1GB for installation, additional space for analysis outputs

### Recommended for GPU Acceleration

- **GPU**: NVIDIA GPU with CUDA support (8GB+ VRAM)
- **CUDA**: Version 11+ with compatible drivers
- **ArrayFire**: Version 3.9.0+ with CUDA backend

## Platform-Specific Notes

### Linux (Ubuntu/Debian)

Install dependencies:
```bash
sudo apt-get update
sudo apt-get install build-essential gfortran cmake \
    libgsl-dev libblas-dev liblapack-dev libopenblas-dev \
    libeigen3-dev libopenmpi-dev openmpi-bin
```

### Linux (CentOS/RHEL)

Install dependencies:
```bash
sudo yum groupinstall "Development Tools"
sudo yum install gsl-devel blas-devel lapack-devel \
    eigen3-devel openmpi-devel
```

### macOS

Install dependencies with Homebrew:
```bash
brew install gsl openblas lapack eigen open-mpi
```

## Troubleshooting

### Common Build Issues

**Missing TEMPO2:**
```bash
# Error: TEMPO2 not found
export TEMPO2=/path/to/tempo2
export PATH=$TEMPO2/bin:$PATH
```

**ArrayFire not found:**
```bash
# Disable GPU support
./configure --without-arrayfire

# Or specify path
./configure --with-arrayfire=/path/to/arrayfire
```

**Eigen3 headers missing:**
```bash
# Auto-detection usually works, but if needed:
./configure --with-eigen=/usr/include/eigen3
```

**MPI compiler issues:**
```bash
# Ensure MPI compilers are in PATH
which mpicc mpicxx

# Or specify explicitly
export CC=mpicc
export CXX=mpicxx
```

### Runtime Issues

**Plugin not found:**
```bash
# Reinstall plugin
make temponest-install

# Check installation path
echo $TEMPO2/plugins/
```

**Segmentation faults:**
```bash
# Use debug build for detailed error information
make clean
./configure --enable-debug
make temponest
```

**JSON parsing errors:**
- Validate JSON syntax with online tools
- Check required fields in configuration
- Ensure file paths are correct

## Environment Setup

### Shell Configuration

Add to your `.bashrc` or `.zshrc`:

```bash
# TempoNest environment
export TEMPO2=/path/to/tempo2
export PATH=$TEMPO2/bin:$PATH

# Optional: GPU acceleration
export ARRAYFIRE_PATH=/path/to/arrayfire
export LD_LIBRARY_PATH=$ARRAYFIRE_PATH/lib64:$LD_LIBRARY_PATH

# Optional: MultiNest
export MULTINEST_PATH=/path/to/multinest
export LD_LIBRARY_PATH=$MULTINEST_PATH:$LD_LIBRARY_PATH
```

### Module Systems (HPC)

For module-based systems:
```bash
# Load required modules
module load tempo2 gsl openblas eigen3 openmpi

# Optional: GPU modules
module load cuda arrayfire
```

## Next Steps

After successful installation:

1. **Configuration**: See [configuration.md](configuration.md) for JSON setup
2. **Examples**: Try examples in [examples.md](examples.md)
3. **Models**: Learn about model elements in [models.md](models.md)
4. **Python Tools**: Explore GUI and plotting tools in [python_tools.md](python_tools.md)

## Support

- **Documentation**: See other files in `docs/` directory
- **Issues**: Report problems with specific error messages
- **Testing**: Use provided test data for verification