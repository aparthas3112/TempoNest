# TempoNest Performance Optimization Action Plan

## Overview
This document outlines practical performance improvements for TempoNest, focusing on GPU optimization, PolyChord integration, and parallel sampling strategies. These improvements target the primary use case of 30-70 parameter sampling where MultiNest becomes inefficient.

## Current Performance Bottlenecks

### Sampling Dimensionality Challenge
- **30 parameters**: ~10⁶ - 10⁷ likelihood evaluations
- **50 parameters**: ~10⁷ - 10⁸ likelihood evaluations  
- **70 parameters**: ~10⁸ - 10⁹ likelihood evaluations

Even with GPU acceleration (1-2ms per evaluation), a 70-parameter problem requires 28 hours to 11 days of computation.

## Implementation Plan

### Phase 1: GPU Implementation Improvements (Week 1)

#### Objectives
- Optimize existing ArrayFire GPU implementation
- Reduce memory allocations and improve parallelism
- Target: 2-3x speedup on current GPU performance

#### Specific Optimizations

1. **Vectorize Diagonal Update**
```cpp
// Current (sequential)
for (int i = 0; i < totCoeff; i++) {
    TNT(start_idx + i, start_idx + i) += pc_inv(i);
}

// Optimized (vectorized)
af::seq diag_idx(start_idx, start_idx + totCoeff - 1);
TNT(diag_idx, diag_idx) += af::diag(pc_inv, 0, false);
```

2. **Persistent Memory Allocation**
- Create class-level cache for frequently used arrays
- Avoid repeated GPU memory allocation/deallocation
- Pre-allocate NT, TNT, and working arrays

3. **Batched Likelihood Evaluation**
- Process multiple parameter sets in single GPU kernel
- Amortize kernel launch overhead
- Optimize for MultiNest's parallel proposals

4. **Remove Unnecessary Synchronization**
- Eliminate `af::sync()` calls except when retrieving final results
- Use asynchronous evaluation where possible

#### Testing Strategy
- Benchmark against current implementation
- Verify numerical accuracy (< 1e-12 difference)
- Test with varying problem sizes (100-5000 observations)

### Phase 2: PolyChord Integration (Weeks 2-3)

#### Why PolyChord?
- **Designed for high dimensions**: Slice sampling scales linearly vs MultiNest's exponential scaling
- **50-parameter efficiency**: PolyChord ~10% vs MultiNest <0.001%
- **Expected speedup**: 10-100x for 30-70 parameters

#### Implementation Details

1. **Use PolyChordLite** (open source version)
   - Core algorithm sufficient for our needs
   - Avoids licensing complications
   - Well-maintained C++/Fortran interface

2. **Sampler Configuration**
```cpp
class polychord_sampler_t : public sampler_t {
    Settings settings;
    
    polychord_sampler_t(int n_dims) {
        settings.nlive = 25 * n_dims;       // Much less than MultiNest
        settings.num_repeats = 5 * n_dims;   // Slice sampling efficiency
        settings.precision_criterion = 0.001;
    }
};
```

3. **Integration Points**
   - Likelihood wrapper (reuse existing GPU implementation)
   - Prior transformations (maintain current format)
   - Output file compatibility
   - JSON configuration support

4. **Configuration Example**
```json
{
  "sampler": {
    "type": "polychord",
    "polychord": {
      "num_live": 1250,      // 25 * 50 parameters
      "num_repeats": 250,    // 5 * 50 parameters
      "precision": 0.001,
      "clustering": true,
      "feedback": 1
    }
  }
}
```

### Phase 3: GPU Stream Parallelism (Week 4)

#### Concept
Run multiple independent sampling chains concurrently on single GPU using different execution streams.

#### Implementation Approach

1. **Multi-Chain MultiNest**
```cpp
class stream_parallel_multinest {
    static constexpr int NUM_STREAMS = 4;
    
    void run_parallel() {
        // Each stream runs independent chain with fewer live points
        // Combine evidence using log-sum-exp
    }
};
```

2. **ArrayFire Strategy**
   - Use thread-level parallelism (OpenMP)
   - Batch operations for multiple chains
   - Automatic work distribution

3. **Expected Benefits**
   - 2.5-3.5x speedup from better GPU utilization
   - More robust evidence estimation
   - Better exploration of multimodal posteriors

## Performance Targets

### Combined Improvements (50-parameter example)
- **Current**: 10-20 hours (MultiNest + GPU)
- **After GPU optimization**: 5-10 hours (2x improvement)
- **With PolyChord**: 30-60 minutes (20x improvement)
- **With stream parallelism**: 15-30 minutes (40x improvement)

### Scaling Expectations
| Parameters | Current Time | Optimized GPU | PolyChord | + Streams |
|------------|-------------|---------------|-----------|-----------|
| 30         | 2-4 hours   | 1-2 hours     | 6-12 min  | 3-6 min   |
| 50         | 10-20 hours | 5-10 hours    | 30-60 min | 15-30 min |
| 70         | 1-2 days    | 12-24 hours   | 2-4 hours | 1-2 hours |

## Success Criteria

1. **GPU Optimization**
   - 2-3x speedup on likelihood evaluation
   - Maintain numerical accuracy (< 1e-12)
   - No increase in memory usage

2. **PolyChord Integration**
   - Successfully sample 50+ parameter models
   - Validate posteriors against MultiNest
   - Evidence agreement within 0.5 log units

3. **Stream Parallelism**
   - 2.5x+ speedup from parallelization
   - Stable evidence combination
   - Efficient GPU utilization (>80%)

## Risk Mitigation

1. **GPU Optimization Risks**
   - Numerical stability: Extensive validation against current implementation
   - Memory limits: Profile memory usage, implement fallbacks

2. **PolyChord Integration Risks**
   - API changes: Pin to specific PolyChord version
   - Prior compatibility: Careful transformation validation

3. **Parallel Sampling Risks**
   - Evidence combination: Use proven statistical methods
   - Random seed management: Ensure independent chains

## Next Steps

1. **Immediate**: Begin GPU optimization implementation
2. **Week 1 checkpoint**: Benchmark GPU improvements
3. **Week 2**: Start PolyChord integration if GPU optimization successful
4. **Week 4**: Implement stream parallelism as final optimization

This plan provides a clear path to 20-40x performance improvement for high-dimensional problems while maintaining scientific accuracy and evidence computation capabilities.