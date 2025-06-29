# TempoNest Documentation

## Overview

TempoNest is a Bayesian pulsar timing analysis tool that integrates with Tempo2 as a plugin. It provides advanced noise modeling, GPU acceleration, and robust parameter estimation for pulsar timing array analyses.

## 🚀 Quick Start

**New to TempoNest?** Start here:

1. **[Installation Guide](installation.md)** - Build and install TempoNest
2. **[Configuration Guide](configuration.md)** - Learn JSON configuration syntax
3. **[Examples](examples.md)** - Working configurations for common analyses
4. **[Python Tools](python_tools.md)** - GUI for configuration and plotting tools

## 📚 Documentation Structure

Our documentation is organized into focused, task-oriented guides:

### 🔧 Getting Started
- **[Installation Guide](installation.md)** - Dependencies, build process, and testing
- **[Configuration Guide](configuration.md)** - JSON syntax, validation, and best practices  
- **[Examples](examples.md)** - Complete working configurations for common scenarios

### 🧮 Scientific Content
- **[Model Elements](models.md)** - All noise models, error models, and physical effects
- **[JUMP Parameters](jumps.md)** - Automatic JUMP handling and marginalization
- **[GPU Optimization](gpu_optimization.md)** - Performance, benchmarks, and ArrayFire details

### 🛠️ Tools and Utilities
- **[Python Tools](python_tools.md)** - Configuration GUI and corner plotting
- **[Troubleshooting](troubleshooting.md)** - Common issues, debugging, and solutions

## 🎯 Use Case Navigation

**I want to...**

### Install and Set Up TempoNest
→ **[Installation Guide](installation.md)**
- System requirements and dependencies
- Build commands and configuration options
- Testing and verification procedures

### Create My First Analysis
→ **[Examples](examples.md)** → Basic Red Noise Analysis
- Copy the basic configuration
- Understand the key parameters
- Run your first analysis

### Understand Model Components
→ **[Model Elements](models.md)**
- Physical noise processes (red noise, DM noise)
- Instrumental effects (EFAC, EQUAD, ECORR) 
- Solar wind modeling
- Mathematical formulations

### Configure Complex Analyses
→ **[Configuration Guide](configuration.md)** + **[Examples](examples.md)**
- JSON syntax and validation
- Parameter types and priors
- Multi-element configurations

### Use the Configuration GUI
→ **[Python Tools](python_tools.md)** → Configuration GUI
- Interactive JSON generation
- Parameter validation and suggestions
- Export and usage instructions

### Analyze and Plot Results
→ **[Python Tools](python_tools.md)** → Corner Plot Tool
- Basic corner plots
- Multi-dataset comparisons
- Parameter mapping and labeling

### Optimize Performance
→ **[GPU Optimization](gpu_optimization.md)**
- Hardware requirements
- Performance benchmarks
- ArrayFire integration details

### Solve Problems
→ **[Troubleshooting](troubleshooting.md)**
- Build and installation issues
- Runtime errors and debugging
- Performance optimization

## 🔬 Scientific Context

### Key Capabilities

**Noise Modeling:**
- Power-law red noise and DM noise
- Comprehensive instrumental effects (EFAC, EQUAD, ECORR)
- Solar wind corrections (deterministic and stochastic)

**Advanced Features:**
- Automatic JUMP parameter marginalization
- GPU acceleration with ArrayFire
- Multi-flag ECORR support
- Legacy compatibility with robust improvements

**Analysis Framework:**
- Nested sampling with MultiNest
- Compositional model system
- Time-span dependent frequency sampling
- Bayesian evidence calculation

### Typical Use Cases

1. **Pulsar Timing Array Analyses** - Individual pulsar noise characterization
2. **Gravitational Wave Studies** - Noise modeling for detection pipelines
3. **Instrumental Systematics** - Multi-backend/frontend correlation studies
4. **Solar Wind Research** - Interplanetary medium effect characterization
5. **Method Development** - Testing new noise models and analysis techniques

## 🏗️ Technical Architecture

**Core Design:**
```
TempoNest Plugin (C++)
├── Model Space (Compositional)
│   ├── Timing Model (Tempo2 integration)
│   ├── Noise Models (Red, DM, Solar wind)
│   └── Error Models (EFAC, EQUAD, ECORR)
├── Likelihood Engine
│   ├── CPU Implementation (Eigen)
│   └── GPU Implementation (ArrayFire)
└── Sampling Framework
    ├── MultiNest Integration
    └── Future: HMC with JAX
```

**Integration:**
- **Tempo2 Plugin**: Seamless integration with existing workflows
- **JSON Configuration**: Modern, validated configuration system
- **Python Tools**: User-friendly interfaces and analysis tools
- **Cross-Platform**: Linux, macOS support with standard dependencies

## 📊 Performance Characteristics

**Scaling:**
- **Parameters**: Efficient for 10-100+ parameter problems
- **GPU Acceleration**: 1-2x speedup vs 300+ CPU cores for large problems
- **Memory**: Proportional to N² (observations + noise coefficients)

**Typical Runtimes:**
- **Testing** (500 live points): 20-60 minutes
- **Standard** (1000 live points): 1-2 hours  
- **Publication** (4000 live points): 4-8 hours

## 🆕 Recent Updates (2025)

### Major Improvements
- **Repository Reorganization**: Modern directory structure with focused documentation
- **GPU Optimizations**: 10-50x speedup for matrix operations through parallel broadcasting
- **Coefficient Calculation**: Scientific time-span dependent frequency sampling
- **ECORR Implementation**: Complete multi-flag epoch correlation support
- **Build System**: Automated build script with comprehensive dependency handling

### Mathematical Verification
- **JUMP Parameters**: Confirmed identical results to legacy marginalization
- **Noise Models**: Verified mathematical consistency across all implementations
- **GPU vs CPU**: Machine-precision agreement for all test cases

## 📖 Documentation Philosophy

Our documentation follows these principles:

1. **Task-Oriented**: Organized by what you want to accomplish
2. **Complete Examples**: Working configurations you can copy and modify
3. **Cross-Referenced**: Links between related topics and concepts
4. **Scientifically Accurate**: Mathematical formulations and physical interpretations
5. **Beginner Friendly**: Clear explanations without assuming prior knowledge
6. **Expert Friendly**: Detailed technical information for advanced users

## 🤝 Getting Help

**Start Here:**
1. Check the relevant documentation section above
2. Look for your specific issue in [Troubleshooting](troubleshooting.md)
3. Try the working examples in [Examples](examples.md)

**For Issues:**
- Include system information, error messages, and reproducible examples
- Check that basic installation and test cases work first
- Provide configuration files and input data when relevant

## 📁 Legacy Documentation

The original comprehensive documentation is preserved in `temponest_doc.md` for reference, but the focused guides above provide better navigation and usability.

---

**Quick Navigation:**
[Install](installation.md) | [Configure](configuration.md) | [Examples](examples.md) | [Models](models.md) | [Tools](python_tools.md) | [GPU](gpu_optimization.md) | [Troubleshoot](troubleshooting.md)