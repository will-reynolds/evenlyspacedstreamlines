#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "enhanced_engine.h"
#include "geometry.h"

//-----------------------------------------------------------------------------
EnhancedStreamlineEngine::EnhancedStreamlineEngine() : StreamlineEngine()
{
    preference_factor = 1.5;
    optimal_mask = nullptr;
    cache_index = 0;
    clear_distance_cache();
}

//-----------------------------------------------------------------------------
EnhancedStreamlineEngine::~EnhancedStreamlineEngine()
{
    delete optimal_mask;
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::setup_enhanced(double radius_, int max_length, int max_nb_seeds_, 
                                              int avoid_u_turns_, double max_angle, 
                                              int oriented_streamlines_,
                                              double singularity_mask_radius,
                                              unsigned int random_seed, int parallel_, int num_threads,
                                              double preference_factor_)
{
    // Call parent setup
    setup(radius_, max_length, max_nb_seeds_, avoid_u_turns_, max_angle, oriented_streamlines_,
          singularity_mask_radius, random_seed, parallel_, num_threads);
    
    preference_factor = preference_factor_;
    
    // Initialize enhanced masking system
    delete optimal_mask;
    optimal_mask = new OptimalSpacingMask(mesh.nt, radius, preference_factor);
    
    // Copy existing binary masks to enhanced system as forbidden regions
    for (int i = 0; i < mesh.nt; i++) {
        for (int j = 0; j < tracer.mask[i].n; j++) {
            double s1 = tracer.mask[i].t[2*j];
            double s2 = tracer.mask[i].t[2*j+1];
            optimal_mask->set_forbidden_mask(i, s1, s2);
        }
    }
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::clear_distance_cache()
{
    for (int i = 0; i < CACHE_SIZE; i++) {
        distance_cache[i].is_valid = false;
    }
    cache_index = 0;
}

//-----------------------------------------------------------------------------
double EnhancedStreamlineEngine::get_cached_distance(int idx, double s)
{
    int hash = (idx * 1000 + (int)(s * 1000)) % CACHE_SIZE;
    
    if (distance_cache[hash].is_valid &&
        distance_cache[hash].triangle_idx == idx &&
        fabs(distance_cache[hash].s_coord - s) < 1e-6) {
        return distance_cache[hash].min_distance;
    }
    
    return -1.0; // Not found
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::cache_distance(int idx, double s, double distance)
{
    int hash = (idx * 1000 + (int)(s * 1000)) % CACHE_SIZE;
    
    distance_cache[hash].triangle_idx = idx;
    distance_cache[hash].s_coord = s;
    distance_cache[hash].min_distance = distance;
    distance_cache[hash].is_valid = true;
}

//-----------------------------------------------------------------------------
double EnhancedStreamlineEngine::compute_min_distance_to_streamlines(int idx, double s)
{
    // Check cache first
    double cached = get_cached_distance(idx, s);
    if (cached >= 0.0) {
        return cached;
    }
    
    double P[3], Q[3];
    mesh.get_segment(idx, s, P, Q);
    
    double min_distance = 1e9; // Large initial value
    
    // Check distance to all existing streamline segments in the spatial data structure
    std::vector<SegmentNode> nearby_segments;
    allsegments.query_neighborhood(idx, nearby_segments);
    
    for (size_t i = 0; i < nearby_segments.size(); i++) {
        const SegmentNode& seg = nearby_segments[i];
        
        // Distance from point P to segment
        double dist_P = point_to_segment_distance(P, seg.P, seg.Q);
        double dist_Q = point_to_segment_distance(Q, seg.P, seg.Q);
        
        // Distance between segments
        double seg_dist = segment_to_segment_distance(P, Q, seg.P, seg.Q);
        
        double current_min = std::min({dist_P, dist_Q, seg_dist});
        
        if (current_min < min_distance) {
            min_distance = current_min;
        }
    }
    
    // If no nearby segments, distance is effectively infinite
    if (nearby_segments.empty()) {
        min_distance = radius * 10.0; // Large distance
    }
    
    // Cache the result
    cache_distance(idx, s, min_distance);
    
    return min_distance;
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::compute_distance_based_mask(int idx_neigh, double* P, double* Q, double radius_)
{
    // Sample multiple points along the triangle's s-coordinate range
    const int num_samples = 20;
    
    for (int i = 0; i <= num_samples; i++) {
        double s = (double)i / (double)num_samples;
        
        double seg_P[3], seg_Q[3];
        mesh.get_segment(idx_neigh, s, seg_P, seg_Q);
        
        // Compute distance from this triangle segment to the new streamline segment PQ
        double dist = segment_to_segment_distance(seg_P, seg_Q, P, Q);
        
        // Add weighted interval based on distance
        double ds = radius_ * mesh.shrink_factor(idx_neigh);
        optimal_mask->set_distance_based_mask(idx_neigh, s - ds, s + ds, dist);
    }
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::update_mask_enhanced(SegmentDeque* queue)
{
    #define QUEUE_FOR_EACH_ELEMENT(queue, idx, s) \
        queue->rewind(); \
        int idx=0; double s=0; \
        while (queue->peek_next(idx, s))

    #define GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q) \
        double P[3], Q[3]; \
        if (queue->is_extremity()) { \
            double t1, t2; \
            queue->get_data(t1, t2); \
            mesh.get_partial_segment(idx, s, t1, t2, P, Q); \
        } else mesh.get_segment(idx, s, P, Q);

    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        // Traditional forbidden masking for the same triangle
        double ds = radius * mesh.shrink_factor(idx);
        optimal_mask->set_forbidden_mask(idx, s - ds, s + ds);
        
        // Enhanced distance-based masking for neighboring triangles
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        int nb_neigh, *neigh_list;
        mesh.get_neighborhood(idx, neigh_list, nb_neigh);
        
        for (int k = 0; k < nb_neigh; k++) {
            int idx_neigh = neigh_list[k];
            compute_distance_based_mask(idx_neigh, P, Q, radius);
        }
        
        // Update the spatial search structure
        allsegments.insert(P, Q, idx);
    }
    
    // Clear distance cache after updating masks
    clear_distance_cache();
}

//-----------------------------------------------------------------------------
void EnhancedStreamlineEngine::update_mask_enhanced_parallel(SegmentDeque* queue)
{
    // Parallel version - similar structure but with OpenMP
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 4)
        for (int pos = queue->first; pos <= queue->last; pos++) {
            int idx = queue->array[pos].idx;
            double s = queue->array[pos].s;
            
            // Same triangle forbidden masking
            double ds = radius * mesh.shrink_factor(idx);
            #pragma omp critical
            {
                optimal_mask->set_forbidden_mask(idx, s - ds, s + ds);
            }
            
            // Distance-based masking for neighbors
            double P[3], Q[3];
            if (pos == queue->first || pos == queue->last) {
                double t1, t2;
                queue->get_data(pos, t1, t2);
                mesh.get_partial_segment(idx, s, t1, t2, P, Q);
            } else {
                mesh.get_segment(idx, s, P, Q);
            }
            
            int nb_neigh, *neigh_list;
            mesh.get_neighborhood(idx, neigh_list, nb_neigh);
            
            for (int k = 0; k < nb_neigh; k++) {
                int idx_neigh = neigh_list[k];
                #pragma omp critical
                {
                    compute_distance_based_mask(idx_neigh, P, Q, radius);
                }
            }
        }
    } // end parallel
    
    // Update spatial search structure (sequential part)
    #define QUEUE_FOR_EACH_ELEMENT(queue, idx, s) \
        queue->rewind(); \
        int idx=0; double s=0; \
        while (queue->peek_next(idx, s))

    #define GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q) \
        double P[3], Q[3]; \
        if (queue->is_extremity()) { \
            double t1, t2; \
            queue->get_data(t1, t2); \
            mesh.get_partial_segment(idx, s, t1, t2, P, Q); \
        } else mesh.get_segment(idx, s, P, Q);

    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        allsegments.insert(P, Q, idx);
    }
    
    clear_distance_cache();
}

//-----------------------------------------------------------------------------
int EnhancedStreamlineEngine::generate_seed_points_enhanced()
{
    int count = optimal_mask->count_admissible_triangles();
    if (count == 0) return 0;
    
    // Generate random numbers
    for (int i = 0; i < 2 * max_nb_seeds; i++)
        rnd_stream[i] = rand();
    
    // Enhanced seed selection with optimal spacing preference
    int nb_seeds = 0;
    std::vector<int> admissible_triangles;
    
    // Build list of triangles that have admissible regions
    for (int i = 0; i < mesh.nt; i++) {
        if (optimal_mask->mask[i].total_weight > 0.0 || optimal_mask->mask[i].is_empty()) {
            admissible_triangles.push_back(i);
        }
    }
    
    if (admissible_triangles.empty()) return 0;
    
    // Select seeds with preference for optimal spacing
    for (int i = 0; i < max_nb_seeds && nb_seeds < max_nb_seeds; i++) {
        int triangle_idx = admissible_triangles[rnd_stream[2*i] % admissible_triangles.size()];
        
        double s = optimal_mask->pick_random_optimal(triangle_idx, rnd_stream[2*i+1]);
        
        if (s >= 0.0) { // Valid seed found
            idx_seed[nb_seeds] = triangle_idx;
            s_seed[nb_seeds] = s;
            nb_seeds++;
        }
    }
    
    return nb_seeds;
}