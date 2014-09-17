#include <armadillo>
#include <cmath>
#include <libwfa/export/grid3d.h>
#include <libwfa/libwfa_exception.h>
#include "grid3d_test.h"

using namespace arma;


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
    if (accu(g.npts()) != 3)
        fail_test(testname, __FILE__, __LINE__, "npts");
    if (fabs(accu(g.origin())) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "origin");
    for (unsigned int i = 0; i < 3; i++) {
        vec ref(3, fill::zeros); ref(i) = 1.;
        if (fabs(accu(ref - g.direction(i))) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "direction");
    }
}


/** \brief Test construction of a symmetric grid
 **/
void grid3d_test::test_1b() {

    static const char *testname = "grid3d_test::test_1b()";

    grid3d g(10, 0.5);

    unsigned int n_[3] = { 10, 10, 10 };
    double x0_[3] = { -2.25, -2.25, -2.25 };

    Col<unsigned int> n_ref(n_, 3);
    vec x0_ref(x0_, 3);
    if (accu(n_ref - g.npts()) != 0)
        fail_test(testname, __FILE__, __LINE__, "npts");
    if (fabs(accu(x0_ref - g.origin())) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "origin");

    for (unsigned int i = 0; i < 3; i++) {
        vec ei_ref(3, fill::zeros);
        ei_ref(i) = 0.5;
        if (fabs(accu(ei_ref - g.direction(i))) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "direction");
    }

}


/** \brief Test construction of a non-symmetric grid
 **/
void grid3d_test::test_1c() {

    static const char *testname = "grid3d_test::test_1c()";

    unsigned int n_[3] = { 10, 5, 8 };
    double d_[3] = { 0.2, 1.0, 0.5 };
    double x0_[3] = { -0.9, -2, -1.75 };

    Col<unsigned int> n(n_, 3);
    vec d(d_, 3), x0_ref(x0_, 3);

    grid3d g(n, d);
    if (accu(n - g.npts()) > 0)
        fail_test(testname, __FILE__, __LINE__, "npts");
    if (fabs(accu(x0_ref - g.origin())) > 1e-15)
        fail_test(testname, __FILE__, __LINE__, "origin");
    for (unsigned int i = 0; i < 3; i++) {
        vec ei(3, fill::zeros);
        ei(i) = d(i);
        if (fabs(accu(ei - g.direction(i))) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "direction");
    }
}


/** \brief Test construction of grid with change of origin
 **/
void grid3d_test::test_2a() {

    static const char *testname = "grid3d_test::test_2a()";

    double x0a_[3] = { -4.5, -4.5, -4.5 };
    double x0b_[3] = { 0., 0., 0. };
    vec x0a_ref(x0a_, 3), x0(x0b_, 3);

    grid3d g(10, 1.0);
    if (fabs(accu(x0a_ref - g.origin())) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "origin");

    g.set_origin(x0);
    if (fabs(accu(x0 - g.origin())) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "origin");
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
    double x0_ref[3] = { 0.0, 0.0, 0.0 };
    for (unsigned int i = 0; i < 3; i++)
        g.set_direction(i, n[i], vec(ev[i], 3));

    if (accu(Col<unsigned int>(n, 3) - g.npts()) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "npts");
    if (fabs(accu(vec(x0_ref, 3) - g.origin())) > 0.0)
        fail_test(testname, __FILE__, __LINE__, "origin");

    for (unsigned int i = 0; i < 3; i++) {
        if (fabs(accu(vec(ev[i], 3) - g.direction(i))) > 0.0)
            fail_test(testname, __FILE__, __LINE__, "direction");
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
    for (unsigned int i = 0; i < 3; i++)
        g.set_direction(i, n[i], vec(ev[i], 3));

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
    double x0_[3] = { 0., 0., 0. };
    vec x0(x0_, 3);
    g.set_origin(x0);

    mat pts(3, 10);
    g.build_pts(23, pts);
    for (unsigned int i = 0, i0 = 0; i < 5; i++) {
        for (unsigned int j = 0; j < 5; j++) {
            for (unsigned int k = 0; k < 5; k++, i0++) {
                if (i0 < 23 || i0 >= 33) continue;

                vec pt_ref = x0 + g.direction(0) * i +
                        g.direction(1) * j + g.direction(2) * k;
                if (fabs(accu(pts.col(i0 - 23) - pt_ref)) > 1e-15)
                    fail_test(testname, __FILE__, __LINE__, "pts");
            }
        }
    }
}

} // namespace libwfa
