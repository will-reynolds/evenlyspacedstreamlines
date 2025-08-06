#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include "preferential_masking.h"

//-----------------------------------------------------------------------------
PreferentialIntervals::PreferentialIntervals(double radius_, double preference_factor_)
{
    radius = radius_;
    preference_factor = preference_factor_;
    n_intervals = 0;
    max_intervals = 8;
    intervals = new WeightedInterval[max_intervals];
    total_weight = 0.0;
}

//-----------------------------------------------------------------------------
PreferentialIntervals::~PreferentialIntervals()
{
    delete [] intervals;
}

//-----------------------------------------------------------------------------
void PreferentialIntervals::clear()
{
    n_intervals = 0;
    total_weight = 0.0;
}

//-----------------------------------------------------------------------------
bool PreferentialIntervals::is_empty()
{
    return (n_intervals == 0);
}

//-----------------------------------------------------------------------------
void PreferentialIntervals::reallocate(int new_max)
{
    WeightedInterval* new_intervals = new WeightedInterval[new_max];
    for (int i = 0; i < n_intervals; i++) {
        new_intervals[i] = intervals[i];
    }
    delete [] intervals;
    intervals = new_intervals;
    max_intervals = new_max;
}

//-----------------------------------------------------------------------------
double PreferentialIntervals::compute_weight(double distance)
{
    if (distance < radius) {
        return 0.0; // Forbidden region
    } else if (distance < radius * preference_factor) {
        // Preferred region: exponentially increasing weight as we get closer to radius
        double normalized_dist = (distance - radius) / (radius * (preference_factor - 1.0));
        return exp(-5.0 * normalized_dist); // Weight decreases exponentially from distance = radius
    } else {
        // Acceptable region: constant low weight
        return 0.1; // Much lower weight than preferred regions
    }
}

//-----------------------------------------------------------------------------
void PreferentialIntervals::add_forbidden_interval(double t1, double t2)
{
    if (t1 >= t2) return;
    
    if (n_intervals >= max_intervals) {
        reallocate(max_intervals * 2);
    }
    
    intervals[n_intervals].t1 = t1;
    intervals[n_intervals].t2 = t2;
    intervals[n_intervals].weight = 0.0;
    intervals[n_intervals].min_distance = 0.0;
    n_intervals++;
    
    merge_intervals();
}

//-----------------------------------------------------------------------------
void PreferentialIntervals::add_weighted_interval(double t1, double t2, double min_distance)
{
    if (t1 >= t2) return;
    
    if (n_intervals >= max_intervals) {
        reallocate(max_intervals * 2);
    }
    
    double weight = compute_weight(min_distance);
    
    intervals[n_intervals].t1 = t1;
    intervals[n_intervals].t2 = t2;
    intervals[n_intervals].weight = weight;
    intervals[n_intervals].min_distance = min_distance;
    n_intervals++;
    
    merge_intervals();
}

//-----------------------------------------------------------------------------
void PreferentialIntervals::merge_intervals()
{
    if (n_intervals <= 1) return;
    
    // Simple merge: combine overlapping intervals with similar weights
    // For more sophisticated implementation, use interval tree
    
    // Sort intervals by t1
    for (int i = 0; i < n_intervals - 1; i++) {
        for (int j = i + 1; j < n_intervals; j++) {
            if (intervals[i].t1 > intervals[j].t1) {
                WeightedInterval temp = intervals[i];
                intervals[i] = intervals[j];
                intervals[j] = temp;
            }
        }
    }
    
    // Update total weight
    total_weight = 0.0;
    for (int i = 0; i < n_intervals; i++) {
        total_weight += intervals[i].weight * (intervals[i].t2 - intervals[i].t1);
    }
}

//-----------------------------------------------------------------------------
double PreferentialIntervals::pick_random_weighted(int rnd)
{
    if (n_intervals == 0) {
        // No constraints, pick uniformly
        return ((double)rnd / (double)RAND_MAX);
    }
    
    if (total_weight <= 0.0) {
        // All intervals are forbidden, return -1
        return -1.0;
    }
    
    // Weighted random selection
    double target = ((double)rnd / (double)RAND_MAX) * total_weight;
    double cumulative_weight = 0.0;
    
    for (int i = 0; i < n_intervals; i++) {
        if (intervals[i].weight <= 0.0) continue;
        
        double interval_weight = intervals[i].weight * (intervals[i].t2 - intervals[i].t1);
        
        if (cumulative_weight + interval_weight >= target) {
            // Select within this interval
            double local_pos = (target - cumulative_weight) / intervals[i].weight;
            double selected_t = intervals[i].t1 + local_pos;
            
            // Clamp to interval bounds
            if (selected_t < intervals[i].t1) selected_t = intervals[i].t1;
            if (selected_t > intervals[i].t2) selected_t = intervals[i].t2;
            
            return selected_t;
        }
        
        cumulative_weight += interval_weight;
    }
    
    // Fallback: select from the last valid interval
    for (int i = n_intervals - 1; i >= 0; i--) {
        if (intervals[i].weight > 0.0) {
            double pos = ((double)(rnd % 1000) / 1000.0);
            return intervals[i].t1 + pos * (intervals[i].t2 - intervals[i].t1);
        }
    }
    
    return -1.0; // No valid intervals
}

//-----------------------------------------------------------------------------
double PreferentialIntervals::get_weight(double t)
{
    for (int i = 0; i < n_intervals; i++) {
        if (t >= intervals[i].t1 && t <= intervals[i].t2) {
            return intervals[i].weight;
        }
    }
    
    // Default weight for unconstrained regions
    return 1.0;
}

//=============================================================================
// OptimalSpacingMask implementation
//=============================================================================

//-----------------------------------------------------------------------------
OptimalSpacingMask::OptimalSpacingMask(int nt_, double radius_, double preference_factor_)
{
    nt = nt_;
    mask = new PreferentialIntervals[nt];
    
    for (int i = 0; i < nt; i++) {
        mask[i] = PreferentialIntervals(radius_, preference_factor_);
    }
}

//-----------------------------------------------------------------------------
OptimalSpacingMask::~OptimalSpacingMask()
{
    delete [] mask;
}

//-----------------------------------------------------------------------------
void OptimalSpacingMask::clear_triangle(int idx)
{
    mask[idx].clear();
}

//-----------------------------------------------------------------------------
void OptimalSpacingMask::set_forbidden_mask(int idx, double s1, double s2)
{
    if (s1 < 0.0) s1 = 0.0;
    if (s2 > 1.0) s2 = 1.0;
    if (s1 < s2) {
        mask[idx].add_forbidden_interval(s1, s2);
    }
}

//-----------------------------------------------------------------------------
void OptimalSpacingMask::set_distance_based_mask(int idx, double s1, double s2, double min_distance)
{
    if (s1 < 0.0) s1 = 0.0;
    if (s2 > 1.0) s2 = 1.0;
    if (s1 < s2) {
        mask[idx].add_weighted_interval(s1, s2, min_distance);
    }
}

//-----------------------------------------------------------------------------
double OptimalSpacingMask::pick_random_optimal(int idx, int rnd)
{
    return mask[idx].pick_random_weighted(rnd);
}

//-----------------------------------------------------------------------------
int OptimalSpacingMask::count_admissible_triangles()
{
    int count = 0;
    for (int i = 0; i < nt; i++) {
        if (mask[i].total_weight > 0.0 || mask[i].is_empty()) {
            count++;
        }
    }
    return count;
}