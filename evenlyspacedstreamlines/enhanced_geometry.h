#ifndef ENHANCED_GEOMETRY_H_
#define ENHANCED_GEOMETRY_H_

// Enhanced geometric distance computations for optimal spacing

// Distance from point to line segment
double point_to_segment_distance(double* point, double* seg_p1, double* seg_p2);

// Distance between two line segments in 3D
double segment_to_segment_distance(double* seg1_p1, double* seg1_p2, 
                                  double* seg2_p1, double* seg2_p2);

// Closest point on segment to given point
void closest_point_on_segment(double* point, double* seg_p1, double* seg_p2, double* closest);

// Parametric distance along segment (0 = seg_p1, 1 = seg_p2)
double parametric_distance_on_segment(double* point, double* seg_p1, double* seg_p2);

#endif