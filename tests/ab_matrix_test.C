#include <libwfa/ab_matrix.h>
#include "ab_matrix_test.h"

namespace libwfa {

using namespace arma;


void ab_matrix_test::perform() throw(libtest::test_exception) {

    test_1a();
    test_1b();
    test_2a();
    test_2b();
    test_3a();
    test_3b();
    test_4a();
    test_4b();
}


void ab_matrix_test::test_1a() {

    static const char *testname = "ab_matrix_test::test_1a()";

    //
    // Test construction of an empty matrix and resizing
    //

    ab_matrix m;
    if (m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta");
    }
    if (m.alpha().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }
    if (m.beta().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "beta != 0");
    }

    Mat<double> &a = m.alpha(), &b = m.beta();
    a.resize(4, 5);
    b.resize(2, 3);
    if (m.nrows_a() != 4 || m.ncols_a() != 5) {
        fail_test(testname, __FILE__, __LINE__, "Size of alpha.");
    }
    if (m.nrows_b() != 2 || m.ncols_b() != 3) {
        fail_test(testname, __FILE__, __LINE__, "Size of beta.");
    }

}


void ab_matrix_test::test_1b() {

    static const char *testname = "ab_matrix_test::test_1b()";

    //
    // Test construction of an empty matrix and resizing
    //

    ab_matrix m(true);
    if (! m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (&(m.alpha()) != &(m.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta");
    }
    if (m.alpha().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }

    m.alpha().resize(3, 4);
    if (m.nrows_a() != m.nrows_b() || m.nrows_a() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# rows.");
    }
    if (m.ncols_a() != m.ncols_b() || m.ncols_a() != 4) {
        fail_test(testname, __FILE__, __LINE__, "# cols.");
    }
}


void ab_matrix_test::test_2a() {

    static const char *testname = "ab_matrix_test::test_2a()";

    //
    // Test construction of non-empty square matrix with alpha == beta
    //

    size_t n = 10;
    ab_matrix m(n);
    if (! m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (m.alpha().n_elem != n * n || m.alpha().n_rows != n) {
        fail_test(testname, __FILE__, __LINE__, "m != 10x10");
    }

    m.beta().reshape(5, 20);
    if (m.alpha().n_elem != 100 || m.alpha().n_rows != 5) {
        fail_test(testname, __FILE__, __LINE__, "m != 5x20");
    }
}


void ab_matrix_test::test_2b() {

    static const char *testname = "ab_matrix_test::test_2b()";

    //
    // Test construction of non-empty matrix with alpha == beta
    //

    ab_matrix m(10, 8);
    if (! m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }

    if (m.alpha().n_elem != 80 || m.alpha().n_rows != 10) {
        fail_test(testname, __FILE__, __LINE__, "m != 10x8");
    }

    m.alpha().resize(4, 10);
    if (m.beta().n_elem != 40 || m.beta().n_rows != 4) {
        fail_test(testname, __FILE__, __LINE__, "m != 4x10");
    }
}


void ab_matrix_test::test_3a() {

    static const char *testname = "ab_matrix_test::test_3a()";

    //
    // Test construction of non-empty matrix with alpha != beta
    //

    ab_matrix m(4, 4, 3);
    if (m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (m.alpha().n_elem != 16 || m.alpha().n_cols != 4) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 4x4");
    }
    if (m.beta().n_elem != 12 || m.beta().n_cols != 4) {
        fail_test(testname, __FILE__, __LINE__, "beta != 3x4");
    }

    m.alpha().resize(2, 3);
    m.beta().reshape(6, 2);
    if (m.alpha().n_elem != 6 || m.alpha().n_cols != 3) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 2x3");
    }
    if (m.beta().n_elem != 12 || m.beta().n_cols != 2) {
        fail_test(testname, __FILE__, __LINE__, "beta != 6x2");
    }
}


void ab_matrix_test::test_3b() {

    static const char *testname = "ab_matrix_test::test_3b()";

    //
    // Test construction of non-empty matrix with alpha != beta
    //

    ab_matrix m(3, 5, 4, 6);
    if (m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (m.alpha().n_elem != 15 || m.alpha().n_cols != 5) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 3x5");
    }
    if (m.beta().n_elem != 24 || m.beta().n_cols != 6) {
        fail_test(testname, __FILE__, __LINE__, "beta != 4x6");
    }

    m.beta().resize(3, 5);
    if (m.alpha().n_elem != 15 || m.alpha().n_cols != 5) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 3x5");
    }
    if (m.beta().n_elem != 15 || m.beta().n_cols != 5) {
        fail_test(testname, __FILE__, __LINE__, "beta != 3x5");
    }
}


void ab_matrix_test::test_4a() {

    static const char *testname = "ab_matrix_test::test_4a()";

    //
    // Test change of non-empty matrix from alpha != beta to alpha == beta
    //

    ab_matrix m(3, 5, 4, 6);
    if (m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }

    // Fill with random numbers
    m.alpha().randu();
    m.beta().randu();

    Mat<double> alpha = m.alpha();
    m.set_alpha_eq_beta();
    if (! m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag).");
    }
    if (&(m.alpha()) != &(m.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }
    if (m.beta().n_elem != 15 || m.beta().n_cols != 5) {
        fail_test(testname, __FILE__, __LINE__, "beta != 3x5");
    }
    if (accu(alpha == m.alpha()) != 15) {
        fail_test(testname, __FILE__, __LINE__, "alpha changed.");
    }
}


void ab_matrix_test::test_4b() {

    static const char *testname = "ab_matrix_test::test_4b()";

    //
    // Test change of non-empty matrix from alpha == beta to alpha != beta
    //

    ab_matrix m(3, 5);
    if (! m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }

    // Fill with random numbers
    m.alpha().randu();
    m.beta().randu();

    m.set_alpha_neq_beta();
    if (m.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta (flag).");
    }
    if (&(m.alpha()) == &(m.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (m.beta().n_elem != 15 || m.beta().n_cols != 5) {
        fail_test(testname, __FILE__, __LINE__, "beta != 3x5");
    }
    if (accu(m.alpha() == m.beta()) != 15) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (contents)");
    }
}


} // namespace libwfa
