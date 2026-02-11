#include "v3math.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_failures = 0;
static const float EPS = 1e-5f;

static void print_v3(const float *v) {
    printf("[%.6f, %.6f, %.6f]", v[0], v[1], v[2]);
}

static void expect_v3(const char *testname, float *actual, float *expected, float tol) {
    if (v3_equals(actual, expected, tol)) {
        printf("PASS: %s\n", testname);
    } else {
        printf("FAIL: %s\n  expected=", testname);
        print_v3(expected);
        printf("\n  actual  =");
        print_v3(actual);
        printf("\n  tol=%.6g\n", tol);
        g_failures++;
    }
}

static void expect_float(const char *testname, float actual, float expected, float tol) {
    if (isnan(expected) && isnan(actual)) {
        printf("PASS: %s\n", testname);
        return;
    }
    if (fabsf(actual - expected) <= tol) {
        printf("PASS: %s\n", testname);
    } else {
        printf("FAIL: %s\n  expected=%.8f actual=%.8f tol=%.6g\n",
               testname, expected, actual, tol);
        g_failures++;
    }
}

// ---- Tests ----

static void test_v3_from_points(void) {
    float a[3] = {1, 2, 3};
    float b[3] = {4, 6, 3};
    float dst[3];

    v3_from_points(dst, a, b);
    float exp1[3] = {3, 4, 0};
    expect_v3("v3_from_points basic", dst, exp1, EPS);

    // overlap: dst == a
    float a2[3] = {1, 2, 3};
    v3_from_points(a2, a2, b);
    expect_v3("v3_from_points overlap dst==a", a2, exp1, EPS);
}

static void test_v3_add(void) {
    float a[3] = {1, -2, 3};
    float b[3] = {4, 5, -6};
    float dst[3];

    v3_add(dst, a, b);
    float exp[3] = {5, 3, -3};
    expect_v3("v3_add basic", dst, exp, EPS);

    // overlap dst==a
    float a2[3] = {1, -2, 3};
    v3_add(a2, a2, b);
    expect_v3("v3_add overlap dst==a", a2, exp, EPS);

    // adding zero vector
    float z[3] = {0,0,0};
    v3_add(dst, a, z);
    expect_v3("v3_add add zero", dst, a, EPS);
}

static void test_v3_subtract(void) {
    float a[3] = {10, 5, -2};
    float b[3] = {3, 7, 4};
    float dst[3];

    v3_subtract(dst, a, b);
    float exp[3] = {7, -2, -6};
    expect_v3("v3_subtract basic", dst, exp, EPS);

    // overlap dst==b (allowed)
    float b2[3] = {3, 7, 4};
    v3_subtract(b2, a, b2);
    expect_v3("v3_subtract overlap dst==b", b2, exp, EPS);

    // subtract from itself -> zero
    float a2[3] = {10, 5, -2};
    v3_subtract(dst, a2, a2);
    float z[3] = {0,0,0};
    expect_v3("v3_subtract self", dst, z, EPS);
}

static void test_v3_dot_product(void) {
    float a[3] = {1, 2, 3};
    float b[3] = {4, -5, 6};
    expect_float("v3_dot_product basic", v3_dot_product(a, b), 12.0f, EPS);

    // orthogonal
    float x[3] = {1,0,0};
    float y[3] = {0,1,0};
    expect_float("v3_dot_product orthogonal", v3_dot_product(x, y), 0.0f, EPS);

    // with itself equals length^2
    expect_float("v3_dot_product self", v3_dot_product(a, a), 14.0f, EPS);
}

static void test_v3_cross_product(void) {
    float x[3] = {1,0,0};
    float y[3] = {0,1,0};
    float dst[3];

    v3_cross_product(dst, x, y);
    float z[3] = {0,0,1};
    expect_v3("v3_cross_product x×y", dst, z, EPS);

    v3_cross_product(dst, y, x);
    float nz[3] = {0,0,-1};
    expect_v3("v3_cross_product y×x", dst, nz, EPS);

    // overlap dst==a
    float a[3] = {1,0,0};
    v3_cross_product(a, a, y);
    expect_v3("v3_cross_product overlap dst==a", a, z, EPS);
}

static void test_v3_scale(void) {
    float v[3] = {1, -2, 3};
    v3_scale(v, 2.0f);
    float exp[3] = {2, -4, 6};
    expect_v3("v3_scale by 2", v, exp, EPS);

    v3_scale(v, 0.5f);
    float exp2[3] = {1, -2, 3};
    expect_v3("v3_scale by 0.5", v, exp2, EPS);

    v3_scale(v, 0.0f);
    float z[3] = {0,0,0};
    expect_v3("v3_scale by 0", v, z, EPS);
}

static void test_v3_length(void) {
    float v[3] = {3,4,12};
    expect_float("v3_length 3-4-12", v3_length(v), 13.0f, 1e-4f);

    float z[3] = {0,0,0};
    expect_float("v3_length zero", v3_length(z), 0.0f, EPS);

    float n[3] = {-1,-2,-2};
    expect_float("v3_length negative components", v3_length(n), 3.0f, 1e-4f);
}

static void test_v3_normalize(void) {
    float v[3] = {3, 0, 4};
    float dst[3];
    v3_normalize(dst, v);
    float exp[3] = {0.6f, 0.0f, 0.8f};
    expect_v3("v3_normalize 3-0-4", dst, exp, 1e-4f);

    // length should be 1
    expect_float("v3_normalize length==1", v3_length(dst), 1.0f, 1e-4f);

    // overlap dst==a
    float v2[3] = {0, 5, 0};
    v3_normalize(v2, v2);
    float exp2[3] = {0, 1, 0};
    expect_v3("v3_normalize overlap dst==a", v2, exp2, 1e-4f);

    // zero vector -> returns zero vector (and prints Error)
    float z[3] = {0,0,0};
    v3_normalize(dst, z);
    float expz[3] = {0,0,0};
    expect_v3("v3_normalize zero vector", dst, expz, EPS);
}

static void test_v3_angle_quick_and_angle(void) {
    float x[3] = {1,0,0};
    float y[3] = {0,1,0};

    // quick returns cos(theta)
    expect_float("v3_angle_quick x,y cos=0", v3_angle_quick(x,y), 0.0f, 1e-5f);
    expect_float("v3_angle x,y pi/2", v3_angle(x,y), (float)(M_PI/2.0), 1e-4f);

    float a[3] = {1,0,0};
    float b[3] = {1,0,0};
    expect_float("v3_angle_quick same cos=1", v3_angle_quick(a,b), 1.0f, 1e-5f);
    expect_float("v3_angle same 0", v3_angle(a,b), 0.0f, 1e-4f);

    float c[3] = {-1,0,0};
    expect_float("v3_angle_quick opposite cos=-1", v3_angle_quick(a,c), -1.0f, 1e-5f);
    expect_float("v3_angle opposite pi", v3_angle(a,c), (float)M_PI, 1e-4f);
}

static void test_v3_reflect(void) {
    // reflect straight down off an upward normal => straight up
    float v[3] = {0, -1, 0};
    float n[3] = {0, 1, 0};
    float dst[3];
    v3_reflect(dst, v, n);
    float exp[3] = {0, 1, 0};
    expect_v3("v3_reflect simple", dst, exp, 1e-5f);

    // non-unit normal should still work
    float n2[3] = {0, 10, 0};
    v3_reflect(dst, v, n2);
    expect_v3("v3_reflect non-unit normal", dst, exp, 1e-5f);

    // overlap dst==v
    float v2[3] = {1, -1, 0};
    float n3[3] = {0, 1, 0};
    v3_reflect(v2, v2, n3);
    float exp2[3] = {1, 1, 0};
    expect_v3("v3_reflect overlap dst==v", v2, exp2, 1e-5f);
}

int main(void) {
    printf("=== v3test: 3D Math Library Unit Tests ===\n\n");

    test_v3_from_points();
    test_v3_add();
    test_v3_subtract();
    test_v3_dot_product();
    test_v3_cross_product();
    test_v3_scale();
    test_v3_length();
    test_v3_normalize();
    test_v3_angle_quick_and_angle();
    test_v3_reflect();

    printf("\n=== Summary ===\n");
    if (g_failures == 0) {
        printf("ALL TESTS PASSED\n");
        return 0;
    } else {
        printf("FAILURES: %d\n", g_failures);
        return 1;
    }
}
