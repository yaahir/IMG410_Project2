#ifndef V3MATH_H
#define V3MATH_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// form v3 from a to b: dst = b - a
void v3_from_points(float *dst, float *a, float *b);

void v3_add(float *dst, float *a, float *b);
void v3_subtract(float *dst, float *a, float *b);

float v3_dot_product(float *a, float *b);

void v3_cross_product(float *dst, float *a, float *b);

// NOTE: signature is given by assignment; this scales the vector in-place:
// dst = dst * s
void v3_scale(float *dst, float s);

float v3_angle(float *a, float *b);        // angle between a and b (radians)
float v3_angle_quick(float *a, float *b);  // returns cos(theta), no acos

// reflect v about n: dst = v - 2*proj_n(v) (n need not be normalized)
void v3_reflect(float *dst, float *v, float *n);

float v3_length(float *a);

void v3_normalize(float *dst, float *a);

// recommended helper for tests
bool v3_equals(float *a, float *b, float tolerance);

#ifdef __cplusplus
}
#endif

#endif
