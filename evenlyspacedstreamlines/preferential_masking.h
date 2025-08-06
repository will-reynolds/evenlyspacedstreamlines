#ifndef PREFERENTIAL_MASKING_H_
#define PREFERENTIAL_MASKING_H_

#include "disjointintervals.h"

//-----------------------------------------------------------------------------
class PreferentialIntervals {
public:
    // Enhanced interval system with distance-based weighting
    // Regions are classified as:
    // - forbidden (distance < radius)
    // - preferred (radius <= distance < radius * preference_factor)  
    // - acceptable (distance >= radius * preference_factor)
    
    struct WeightedInterval {
        double t1, t2;      // interval bounds
        double weight;      // selection weight (0 = forbidden, higher = preferred)
        double min_distance; // minimum distance to existing streamlines
    };
    
    int n_intervals;
    WeightedInterval* intervals;
    double total_weight;
    double radius;
    double preference_factor; // typically 1.5-2.0
    
    PreferentialIntervals(double radius_, double preference_factor_ = 1.5);
    ~PreferentialIntervals();
    
    void clear();
    bool is_empty();
    
    // Add forbidden region (weight = 0)
    void add_forbidden_interval(double t1, double t2);
    
    // Add weighted region based on minimum distance to existing streamlines
    void add_weighted_interval(double t1, double t2, double min_distance);
    
    // Pick random point with distance-based weighting
    double pick_random_weighted(int rnd);
    
    // Get weight at specific coordinate (for debugging)
    double get_weight(double t);
    
private:
    int max_intervals;
    void reallocate(int new_max);
    void merge_intervals();
    double compute_weight(double distance);
};

//-----------------------------------------------------------------------------
// Enhanced triangle masking with distance-based weights
class OptimalSpacingMask {
public:
    PreferentialIntervals* mask;
    int nt; // number of triangles
    
    OptimalSpacingMask(int nt_, double radius_, double preference_factor_ = 1.5);
    ~OptimalSpacingMask();
    
    void clear_triangle(int idx);
    void set_forbidden_mask(int idx, double s1, double s2);
    void set_distance_based_mask(int idx, double s1, double s2, double min_distance);
    
    double pick_random_optimal(int idx, int rnd);
    int count_admissible_triangles();
};

#endif