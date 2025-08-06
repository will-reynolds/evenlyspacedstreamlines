#ifndef ENHANCED_ENGINE_H_
#define ENHANCED_ENGINE_H_

#include "engine.h"
#include "preferential_masking.h"
#include "segmentsearch.h"

//-----------------------------------------------------------------------------
class EnhancedStreamlineEngine : public StreamlineEngine {
public:
    double preference_factor; // Spacing preference factor (1.5-2.0)
    OptimalSpacingMask* optimal_mask;
    
    EnhancedStreamlineEngine();
    ~EnhancedStreamlineEngine();

    // Override setup to include preference factor
    void setup_enhanced(double radius_, int max_length, int max_nb_seeds_, 
                       int avoid_u_turns, double max_angle, 
                       int oriented_streamlines_,
                       double singularity_mask_radius,
                       unsigned int random_seed, int parallel_, int num_threads,
                       double preference_factor_ = 1.5);
    
    // Enhanced masking functions
    void update_mask_enhanced(SegmentDeque* queue);
    void update_mask_enhanced_parallel(SegmentDeque* queue);
    
    // Enhanced seed generation with optimal spacing preference
    int generate_seed_points_enhanced();
    
    // Distance calculation helpers
    double compute_min_distance_to_streamlines(int idx, double s);
    void compute_distance_based_mask(int idx_neigh, double* P, double* Q, double radius);

private:
    // Distance computation cache for efficiency
    struct DistanceCache {
        int triangle_idx;
        double s_coord;
        double min_distance;
        bool is_valid;
    };
    
    static const int CACHE_SIZE = 1024;
    DistanceCache distance_cache[CACHE_SIZE];
    int cache_index;
    
    void clear_distance_cache();
    double get_cached_distance(int idx, double s);
    void cache_distance(int idx, double s, double distance);
};

#endif