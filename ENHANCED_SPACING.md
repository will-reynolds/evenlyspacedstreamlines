# Enhanced Spacing for Evenly-Spaced Streamlines

This document describes the enhanced masking logic that encourages spacing to be above but as close to the radius as possible.

## Overview

The original algorithm uses binary masking:
- **Forbidden regions**: Distance < radius (weight = 0)
- **Allowed regions**: Distance ≥ radius (equal probability)

The enhanced algorithm uses **graduated masking** with distance-based preferences:
- **Forbidden regions**: Distance < radius (weight = 0)
- **Preferred regions**: radius ≤ distance < radius × preference_factor (high weight)
- **Acceptable regions**: Distance ≥ radius × preference_factor (low weight)

## Key Components

### 1. PreferentialIntervals Class

Replaces binary `DisjointIntervals` with weighted interval selection:

```cpp
class PreferentialIntervals {
    struct WeightedInterval {
        double t1, t2;          // Interval bounds
        double weight;          // Selection probability weight
        double min_distance;    // Distance to nearest streamlines
    };
    
    double pick_random_weighted(int rnd);  // Distance-based selection
};
```

**Weight Function:**
```cpp
double compute_weight(double distance) {
    if (distance < radius) {
        return 0.0;  // Forbidden
    } else if (distance < radius * preference_factor) {
        double normalized = (distance - radius) / (radius * (preference_factor - 1.0));
        return exp(-5.0 * normalized);  // Exponentially decreasing from radius
    } else {
        return 0.1;  // Low constant weight for far regions
    }
}
```

### 2. Enhanced Distance Computation

The enhanced system computes actual distances to existing streamlines:

```cpp
double compute_min_distance_to_streamlines(int idx, double s) {
    // Get 3D coordinates of candidate location
    mesh.get_segment(idx, s, P, Q);
    
    // Query spatial data structure for nearby streamlines
    allsegments.query_neighborhood(idx, nearby_segments);
    
    // Compute minimum distance to all nearby segments
    double min_distance = segment_to_segment_distance(P, Q, seg.P, seg.Q);
}
```

### 3. Enhanced Masking Pipeline

The updated masking process:

```cpp
void update_mask_enhanced(SegmentDeque* queue) {
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        // 1. Forbidden masking (same triangle)
        optimal_mask->set_forbidden_mask(idx, s - ds, s + ds);
        
        // 2. Distance-based masking (neighboring triangles)
        for (int k = 0; k < nb_neigh; k++) {
            compute_distance_based_mask(idx_neigh, P, Q, radius);
        }
        
        // 3. Update spatial search structure
        allsegments.insert(P, Q, idx);
    }
}
```

## Usage

### C++ Integration

1. Include the enhanced headers:
```cpp
#include "enhanced_engine.h"
#include "preferential_masking.h"
```

2. Use the enhanced engine:
```cpp
EnhancedStreamlineEngine engine;
engine.setup_enhanced(radius, max_length, max_nb_seeds, 
                     avoid_u_turns, max_angle, oriented_streamlines,
                     singularity_mask_radius, random_seed, parallel, num_threads,
                     preference_factor);  // New parameter

engine.run();  // Uses enhanced masking internally
```

### Python Interface

```python
from evenlyspacedstreamlines import optimal_spaced_streamlines

streamlines, info = optimal_spaced_streamlines(
    vertices, triangles, orientations, radius,
    preference_factor=1.8  # Strong preference for closer spacing
)

print(f"Spacing efficiency: {info['spacing_efficiency']:.3f}")
print(f"Mean spacing: {info['mean_spacing']:.3f}")
```

## Parameter Tuning

### preference_factor

- **1.0**: No preference (equivalent to original algorithm)
- **1.5**: Moderate preference for closer spacing
- **2.0**: Strong preference for closer spacing  
- **>2.5**: Very aggressive spacing (may reduce coverage)

### Weight Function Tuning

The exponential weight function can be adjusted:

```cpp
// More aggressive preference (steeper falloff)
return exp(-10.0 * normalized);

// Gentler preference (gradual falloff)  
return exp(-2.0 * normalized);

// Linear falloff
return 1.0 - normalized;
```

## Performance Considerations

### Computational Overhead

- **Distance computations**: O(k) per candidate location, where k = nearby streamlines
- **Weighted sampling**: O(n) per selection, where n = number of weighted intervals
- **Caching**: Distance computations are cached to reduce redundancy

### Memory Usage

- **Weighted intervals**: ~5x memory vs binary intervals
- **Distance cache**: ~4KB for typical cache size
- **Spatial queries**: Dependent on streamline density

### Optimization Strategies

1. **Adaptive sampling**: Reduce distance computation frequency in dense regions
2. **Hierarchical caching**: Cache at multiple spatial scales
3. **Early termination**: Stop distance queries when minimum weight is reached

## Expected Results

The enhanced algorithm produces:

1. **Improved spacing uniformity**: Mean spacing closer to target radius
2. **Reduced over-dense regions**: Fewer areas with spacing >> radius  
3. **Better coverage efficiency**: More streamlines within the same total area
4. **Controlled trade-offs**: Tunable preference strength via preference_factor

## Mathematical Foundation

The spacing preference is based on the **probability density function**:

```
P(location) ∝ w(d(location, existing_streamlines))
```

Where:
- `d(·,·)` = minimum distance function
- `w(·)` = weight function (exponential decay from radius)

This creates a **spatially adaptive** sampling distribution that naturally encourages optimal spacing.

## Integration with Existing Code

The enhanced system is designed as a **drop-in replacement**:

1. **Binary compatibility**: Can fall back to original algorithm (preference_factor = 1.0)
2. **API compatibility**: Same function signatures with optional enhancement parameters  
3. **Performance scaling**: Overhead is proportional to preference_factor strength

## Future Enhancements

1. **Anisotropic preferences**: Different spacing preferences along vs across the flow
2. **Adaptive preference factors**: Automatically tune based on local density
3. **Multi-scale optimization**: Hierarchical spacing optimization
4. **Quality metrics**: Automatic evaluation of spacing uniformity