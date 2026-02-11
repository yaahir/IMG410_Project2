#include "v3math.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// ---------- internal helpers ----------
static void v3_error(const char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
}

static bool v3_valid_ptr3(float *p) {
    return p != NULL;
}

static float clampf(float x, float lo, float hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

// ---------- required functions ----------
void v3_from_points(float *dst, float *a, float *b) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_from_points received NULL pointer");
        return;
    }
    // overlap-safe
    float ax = a[0], ay = a[1], az = a[2];
    float bx = b[0], by = b[1], bz = b[2];
    dst[0] = bx - ax;
    dst[1] = by - ay;
    dst[2] = bz - az;
}

void v3_add(float *dst, float *a, float *b) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_add received NULL pointer");
        return;
    }
    float ax = a[0], ay = a[1], az = a[2];
    float bx = b[0], by = b[1], bz = b[2];
    dst[0] = ax + bx;
    dst[1] = ay + by;
    dst[2] = az + bz;
}

void v3_subtract(float *dst, float *a, float *b) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_subtract received NULL pointer");
        return;
    }
    float ax = a[0], ay = a[1], az = a[2];
    float bx = b[0], by = b[1], bz = b[2];
    dst[0] = ax - bx;
    dst[1] = ay - by;
    dst[2] = az - bz;
}

float v3_dot_product(float *a, float *b) {
    if (!v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_dot_product received NULL pointer");
        return NAN;
    }
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void v3_cross_product(float *dst, float *a, float *b) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_cross_product received NULL pointer");
        return;
    }
    // overlap-safe: store inputs first
    float ax = a[0], ay = a[1], az = a[2];
    float bx = b[0], by = b[1], bz = b[2];

    dst[0] = ay*bz - az*by;
    dst[1] = az*bx - ax*bz;
    dst[2] = ax*by - ay*bx;
}

void v3_scale(float *dst, float s) {
    if (!v3_valid_ptr3(dst)) {
        v3_error("v3_scale received NULL pointer");
        return;
    }
    dst[0] *= s;
    dst[1] *= s;
    dst[2] *= s;
}

float v3_length(float *a) {
    if (!v3_valid_ptr3(a)) {
        v3_error("v3_length received NULL pointer");
        return NAN;
    }
    // hypotf for stability
    float x = a[0], y = a[1], z = a[2];
    return sqrtf(x*x + y*y + z*z);
}

void v3_normalize(float *dst, float *a) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(a)) {
        v3_error("v3_normalize received NULL pointer");
        return;
    }
    float x = a[0], y = a[1], z = a[2];
    float len = sqrtf(x*x + y*y + z*z);

    if (len == 0.0f || !isfinite(len)) {
        v3_error("v3_normalize cannot normalize zero-length or non-finite vector");
        // safest behavior: output zero vector
        dst[0] = 0.0f; dst[1] = 0.0f; dst[2] = 0.0f;
        return;
    }

    float inv = 1.0f / len;
    dst[0] = x * inv;
    dst[1] = y * inv;
    dst[2] = z * inv;
}

float v3_angle_quick(float *a, float *b) {
    if (!v3_valid_ptr3(a) || !v3_valid_ptr3(b)) {
        v3_error("v3_angle_quick received NULL pointer");
        return NAN;
    }
    float la = v3_length(a);
    float lb = v3_length(b);
    if (la == 0.0f || lb == 0.0f || !isfinite(la) || !isfinite(lb)) {
        v3_error("v3_angle_quick undefined for zero-length or non-finite vectors");
        return NAN;
    }
    float cosv = v3_dot_product(a, b) / (la * lb);
    // clamp to avoid slight floating error pushing outside [-1, 1]
    return clampf(cosv, -1.0f, 1.0f);
}

float v3_angle(float *a, float *b) {
    float cosv = v3_angle_quick(a, b);
    if (!isfinite(cosv)) return NAN;
    return acosf(cosv);
}

void v3_reflect(float *dst, float *v, float *n) {
    if (!v3_valid_ptr3(dst) || !v3_valid_ptr3(v) || !v3_valid_ptr3(n)) {
        v3_error("v3_reflect received NULL pointer");
        return;
    }

    // Normalize n (but do it overlap-safe)
    float nn[3];
    v3_normalize(nn, n);

    // If normalization failed, nn is zero vector; reflection is undefined.
    float nn_len = v3_length(nn);
    if (nn_len == 0.0f || !isfinite(nn_len)) {
        v3_error("v3_reflect undefined because normal vector is zero-length");
        // reasonable fallback: copy v into dst
        dst[0] = v[0]; dst[1] = v[1]; dst[2] = v[2];
        return;
    }

    float vx = v[0], vy = v[1], vz = v[2];
    float dotvn = vx*nn[0] + vy*nn[1] + vz*nn[2];

    // r = v - 2*dot(v,n)*n
    dst[0] = vx - 2.0f * dotvn * nn[0];
    dst[1] = vy - 2.0f * dotvn * nn[1];
    dst[2] = vz - 2.0f * dotvn * nn[2];
}

bool v3_equals(float *a, float *b, float tolerance) {
    if (!v3_valid_ptr3(a) || !v3_valid_ptr3(b)) return false;
    if (tolerance < 0.0f) tolerance = -tolerance;

    return (fabsf(a[0] - b[0]) <= tolerance) &&
           (fabsf(a[1] - b[1]) <= tolerance) &&
           (fabsf(a[2] - b[2]) <= tolerance);
}
