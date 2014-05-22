#include <libwfa/core/ab_vector.h>
#include "ab_vector_test.h"

namespace libwfa {

using namespace arma;


void ab_vector_test::perform() throw(libtest::test_exception) {

    test_1a();
    test_1b();
    test_2();
    test_3();
    test_4a();
    test_4b();
}


void ab_vector_test::test_1a() {

    static const char *testname = "ab_vector_test::test_1a()";

    //
    // Test construction of an empty vector and resizing
    //

    ab_vector v;
    if (v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta");
    }
    if (v.alpha().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }
    if (v.beta().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "beta != 0");
    }

    Col<double> &a = v.alpha(), &b = v.beta();
    a.resize(4);
    b.resize(2);
    if (v.nrows_a() != 4) {
        fail_test(testname, __FILE__, __LINE__, "Size of alpha.");
    }
    if (v.nrows_b() != 2) {
        fail_test(testname, __FILE__, __LINE__, "Size of beta.");
    }

}


void ab_vector_test::test_1b() {

    static const char *testname = "ab_vector_test::test_1b()";

    //
    // Test construction of an empty vector and resizing
    //

    ab_vector v(true);
    if (! v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (&(v.alpha()) != &(v.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta");
    }
    if (v.alpha().n_elem != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }

    v.alpha().resize(3);
    if (v.nrows_a() != v.nrows_b() || v.nrows_a() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# rows.");
    }
}


void ab_vector_test::test_2() {

    static const char *testname = "ab_vector_test::test_2()";

    //
    // Test construction of non-empty vector with alpha == beta
    //

    size_t n = 10;
    ab_vector v(n);
    if (! v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (v.alpha().n_rows != n) {
        fail_test(testname, __FILE__, __LINE__, "# rows != 10");
    }

    v.beta().resize(5);
    if (v.alpha().n_rows != 5) {
        fail_test(testname, __FILE__, __LINE__, "# rows != 5");
    }
}


void ab_vector_test::test_3() {

    static const char *testname = "ab_vector_test::test_3()";

    //
    // Test construction of non-empty matrix with alpha != beta
    //

    ab_vector v(4, 3);
    if (v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (v.alpha().n_rows != 4) {
        fail_test(testname, __FILE__, __LINE__, "# alpha rows != 4.");
    }
    if (v.beta().n_rows != 3) {
        fail_test(testname, __FILE__, __LINE__, "# beta rows != 3.");
    }

    v.alpha().resize(3);
    if (v.alpha().n_rows != 3) {
        fail_test(testname, __FILE__, __LINE__, "# alpha rows != 3.");
    }
}


void ab_vector_test::test_4a() {

    static const char *testname = "ab_vector_test::test_4a()";

    //
    // Test change of non-empty matrix from alpha != beta to alpha == beta
    //

    ab_vector v(3, 5);
    if (v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }

    // Fill with random numbers
    v.alpha().randu();
    v.beta().randu();

    Col<double> alpha = v.alpha();
    v.set_alpha_eq_beta();
    if (! v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag).");
    }
    if (&(v.alpha()) != &(v.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }
    if (v.beta().n_rows != 3) {
        fail_test(testname, __FILE__, __LINE__, "# beta rows != 3");
    }
    if (accu(alpha == v.alpha()) != v.alpha().n_elem) {
        fail_test(testname, __FILE__, __LINE__, "alpha changed.");
    }
}


void ab_vector_test::test_4b() {

    static const char *testname = "ab_vector_test::test_4b()";

    //
    // Test change of non-empty matrix from alpha == beta to alpha != beta
    //

    size_t n = 5;
    ab_vector v(n);
    if (! v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }

    // Fill with random numbers
    v.alpha().randu();
    v.beta().randu();

    v.set_alpha_neq_beta();
    if (v.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta (flag).");
    }
    if (&(v.alpha()) == &(v.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (v.beta().n_rows != n) {
        fail_test(testname, __FILE__, __LINE__, "# beta rows != 5");
    }
    if (accu(v.alpha() == v.beta()) != v.alpha().n_elem) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (contents)");
    }
}


} // namespace libwfa
