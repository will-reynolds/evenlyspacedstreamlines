#include <math.h>
#include <algorithm>
#include "enhanced_geometry.h"
#include "algebra.h" // For vector operations

//-----------------------------------------------------------------------------
double point_to_segment_distance(double* point, double* seg_p1, double* seg_p2)
{
    double segment[3], point_vec[3];
    vdiff(segment, seg_p2, seg_p1);      // segment = p2 - p1
    vdiff(point_vec, point, seg_p1);     // point_vec = point - p1
    
    double seg_length_sq = vnorm2(segment);
    
    if (seg_length_sq < 1e-12) {
        // Degenerate segment, just distance to point
        return vdist(point, seg_p1);
    }
    
    double t = vdot(point_vec, segment) / seg_length_sq;
    
    // Clamp t to [0, 1] to stay on the segment
    t = std::max(0.0, std::min(1.0, t));
    
    // Closest point on segment
    double closest[3];
    vcopyadd(closest, seg_p1, t, segment);
    
    return vdist(point, closest);
}

//-----------------------------------------------------------------------------
void closest_point_on_segment(double* point, double* seg_p1, double* seg_p2, double* closest)
{
    double segment[3], point_vec[3];
    vdiff(segment, seg_p2, seg_p1);
    vdiff(point_vec, point, seg_p1);
    
    double seg_length_sq = vnorm2(segment);
    
    if (seg_length_sq < 1e-12) {
        // Degenerate segment
        vcopy(closest, seg_p1);
        return;
    }
    
    double t = vdot(point_vec, segment) / seg_length_sq;
    t = std::max(0.0, std::min(1.0, t));
    
    vcopyadd(closest, seg_p1, t, segment);
}

//-----------------------------------------------------------------------------
double parametric_distance_on_segment(double* point, double* seg_p1, double* seg_p2)
{
    double segment[3], point_vec[3];
    vdiff(segment, seg_p2, seg_p1);
    vdiff(point_vec, point, seg_p1);
    
    double seg_length_sq = vnorm2(segment);
    
    if (seg_length_sq < 1e-12) {
        return 0.0; // Degenerate segment
    }
    
    double t = vdot(point_vec, segment) / seg_length_sq;
    return std::max(0.0, std::min(1.0, t));
}

//-----------------------------------------------------------------------------
double segment_to_segment_distance(double* seg1_p1, double* seg1_p2, 
                                  double* seg2_p1, double* seg2_p2)
{
    // Implementation of segment-to-segment distance in 3D
    // Based on the algorithm from "Real-Time Collision Detection" by Christer Ericson
    
    double d1[3], d2[3], r[3];
    vdiff(d1, seg1_p2, seg1_p1);  // Direction of segment 1
    vdiff(d2, seg2_p2, seg2_p1);  // Direction of segment 2  
    vdiff(r, seg1_p1, seg2_p1);   // Vector between start points
    
    double a = vdot(d1, d1);       // Squared length of segment 1
    double e = vdot(d2, d2);       // Squared length of segment 2
    double f = vdot(d2, r);
    
    // Check for degenerate cases
    if (a <= 1e-12 && e <= 1e-12) {
        // Both segments are points
        return vdist(seg1_p1, seg2_p1);
    }
    
    double s, t;
    
    if (a <= 1e-12) {
        // Segment 1 is a point
        s = 0.0;
        t = f / e;
        t = std::max(0.0, std::min(1.0, t));
    } else {
        double c = vdot(d1, r);
        
        if (e <= 1e-12) {
            // Segment 2 is a point
            t = 0.0;
            s = -c / a;
            s = std::max(0.0, std::min(1.0, s));
        } else {
            // General case
            double b = vdot(d1, d2);
            double denom = a * e - b * b;
            
            if (denom != 0.0) {
                s = (b * f - c * e) / denom;
                s = std::max(0.0, std::min(1.0, s));
            } else {
                s = 0.0; // Segments are parallel
            }
            
            t = (b * s + f) / e;
            
            if (t < 0.0) {
                t = 0.0;
                s = -c / a;
                s = std::max(0.0, std::min(1.0, s));
            } else if (t > 1.0) {
                t = 1.0;
                s = (b - c) / a;
                s = std::max(0.0, std::min(1.0, s));
            }
        }
    }
    
    // Compute closest points
    double closest1[3], closest2[3];
    vcopyadd(closest1, seg1_p1, s, d1);
    vcopyadd(closest2, seg2_p1, t, d2);
    
    return vdist(closest1, closest2);
}