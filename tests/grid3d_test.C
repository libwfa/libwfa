#include <cmath>
#include <libwfa/grid3d.h>
#include <libwfa/libwfa_exception.h>
#include "grid3d_test.h"

namespace libwfa {

void grid3d_test::perform() throw(libtest::test_exception) {

    test_1a();
    test_1b();
    test_1c();
    test_2a();
    test_2b();
    test_3();
    test_4();
}


/** \brief Test construction of an empty grid
 **/
void grid3d_test::test_1a() {

    static const char *testname = "grid3d_test::test_1a()";

    grid3d g;
    for (unsigned int i = 0; i < 3; i++) {
        if (g.npts(i) != 1) fail_test(testname, __FILE__, __LINE__, "npts");
        if (fabs(g.origin(i)) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "origin");
        for (unsigned int j = 0; j < 3; j++) {
            double ref = (i == j ? 1. : 0.0);
            if (fabs(g.direction(i, j) - ref) > 0.0)
                fail_test(testname, __FILE__, __LINE__, "direction");
        }
    }
}


/** \brief Test construction of a symmetric grid
 **/
void grid3d_test::test_1b() {

    static const char *testname = "grid3d_test::test_1b()";

    grid3d g(10, 0.5);
    for (unsigned int i = 0; i < 3; i++) {
        if (g.npts(i) != 10) fail_test(testname, __FILE__, __LINE__, "npts");
        if (fabs(g.origin(i) + 2.5) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "origin");
        for (unsigned int j = 0; j < 3; j++) {
            double ref = (i == j ? 0.5 : 0.0);
            if (fabs(g.direction(i, j) - ref) > 0.0)
                fail_test(testname, __FILE__, __LINE__, "direction");
        }
    }
}


/** \brief Test construction of a non-symmetric grid
 **/
void grid3d_test::test_1c() {

    static const char *testname = "grid3d_test::test_1c()";

    unsigned int n[3] = { 10, 5, 8 };
    double d[3] = { 0.2, 1.0, 0.5 };
    grid3d g(n, d);
    for (unsigned int i = 0; i < 3; i++) {
        if (g.npts(i) != n[i]) fail_test(testname, __FILE__, __LINE__, "npts");
        if (fabs(g.origin(i) + d[i] * n[i] / 2.) > 1e-15)
            fail_test(testname, __FILE__, __LINE__, "origin");
        for (unsigned int j = 0; j < 3; j++) {
            double ref = (i == j ? d[i] : 0.0);
            if (fabs(g.direction(i, j) - ref) > 0.0)
                fail_test(testname, __FILE__, __LINE__, "direction");
        }
    }
}


/** \brief Test construction of grid with change of origin
 **/
void grid3d_test::test_2a() {

    static const char *testname = "grid3d_test::test_2a()";

    double x0[3] = { 0., 0., 0. };
    grid3d g(10, 1.0);
    for (unsigned int i = 0; i < 3; i++) {
        if (fabs(g.origin(i) - 5.0) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "origin");
    }

    g.set_origin(x0);
    for (unsigned int i = 0; i < 3; i++) {
        if (fabs(g.origin(i) - x0[i]) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "origin");
    }
}


/** \brief Test construction of grid with new directions
 **/
void grid3d_test::test_2b() {

    static const char *testname = "grid3d_test::test_2b()";


    grid3d g;
    unsigned int n[3] = { 20, 10, 16 };
    double ev[3][3] = {
            { 0.0, 0.2, 0.0 },
            { 0.0, 0.0, 0.5 },
            { 0.2, 0.1, 0.0 } };
    for (unsigned int i = 0; i < 3; i++) g.set_direction(i, n[i], ev[i]);

    for (unsigned int i = 0; i < 3; i++) {
        if (g.npts(i) != n[i]) fail_test(testname, __FILE__, __LINE__, "npts");
        if (fabs(g.origin(i)) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "origin");
        for (unsigned int j = 0; j < 3; j++) {
            if (fabs(g.direction(i, j) - ev[i][j]) > 0.0)
                fail_test(testname, __FILE__, __LINE__, "direction");
        }
    }
}


/** \brief Test check of illegal grid
 **/
void grid3d_test::test_3() {

    static const char *testname = "grid3d_test::test_3()";

#ifdef LIBWFA_DEBUG
    grid3d g;
    unsigned int n[3] = { 20, 10, 16 };
    double ev[3][3] = {
            { 0.0, 0.1, 0.0 },
            { 0.0, 0.0, 0.5 },
            { 0.0, 0.2, 0.0 } };
    for (unsigned int i = 0; i < 3; i++) g.set_direction(i, n[i], ev[i]);

    bool fail = false;
    try {
        g.check();
    } catch (std::exception &e) {
        fail = true;
    }

    if (! fail) {
        fail_test(testname, __FILE__, __LINE__, "check");
    }
#endif
}


/** \brief Test export of grid points
 **/
void grid3d_test::test_4() {

    static const char *testname = "grid3d_test::test_4()";

    grid3d g(5, 0.1);
    double x0[3] = { 0., 0., 0. };
    g.set_origin(x0);

    double pts[30];
    g.build_pts(pts, 0, 5);
    for (unsigned int i = 0, pos = 0; i < 2; i++) {
        for (unsigned int j = 0; j < 5; j++, pos += 3) {
            if (fabs(pts[pos] - x0[0]) > 1e-15)
                fail_test(testname, __FILE__, __LINE__, "pts(0)");
            if (fabs(pts[pos + 1] - (x0[0] + i * 0.1)) > 1e-15)
                fail_test(testname, __FILE__, __LINE__, "pts(1)");
            if (fabs(pts[pos + 2] - (x0[0] + i * 0.1 + j * 0.1)) > 1e-15)
                fail_test(testname, __FILE__, __LINE__, "pts(2)");
        }
    }
}


} // namespace libwfa
