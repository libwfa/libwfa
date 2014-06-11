#include <libwfa/core/ab_orbital_selector.h>
#include "ab_orbital_selector_test.h"

namespace libwfa {


void ab_orbital_selector_test::perform() throw(libtest::test_exception) {

    test_1a();
    test_1b();
    test_2();
    test_3();
    test_4a();
    test_4b();
}


void ab_orbital_selector_test::test_1a() {

    static const char *testname = "ab_orbital_selector_test::test_1a()";

    //
    // Test construction of an empty selector
    //

    ab_orbital_selector s;
    if (s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta");
    }
    if (s.alpha().n_indexes() != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }
    if (s.beta().n_indexes() != 0) {
        fail_test(testname, __FILE__, __LINE__, "beta != 0");
    }

    orbital_selector &sa = s.alpha(), &sb = s.beta();
    sa = orbital_selector(4);
    sb = orbital_selector(3);
    if (s.nidx_a() != 4) {
        fail_test(testname, __FILE__, __LINE__, "Size of alpha.");
    }
    if (s.nidx_b() != 3) {
        fail_test(testname, __FILE__, __LINE__, "Size of beta.");
    }

}


void ab_orbital_selector_test::test_1b() {

    static const char *testname = "ab_orbital_selector_test::test_1b()";

    //
    // Test construction of an empty selector and resizing
    //

    ab_orbital_selector s(true);
    if (! s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (&(s.alpha()) != &(s.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta");
    }
    if (s.alpha().n_indexes() != 0) {
        fail_test(testname, __FILE__, __LINE__, "alpha != 0");
    }

    s.alpha() = orbital_selector(3);
    if (s.nidx_a() != s.nidx_b() || s.nidx_a() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# indexes.");
    }
}


void ab_orbital_selector_test::test_2() {

    static const char *testname = "ab_orbital_selector_test::test_2()";

    //
    // Test construction of non-empty selector with alpha == beta
    //

    size_t n = 10;

    ab_orbital_selector s(n);
    if (! s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag)");
    }
    if (s.alpha().n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes != 10");
    }

    s.beta() = orbital_selector(5);
    if (s.alpha().n_indexes() != 5) {
        fail_test(testname, __FILE__, __LINE__, "# indexes != 5");
    }
}


void ab_orbital_selector_test::test_3() {

    static const char *testname = "ab_orbital_selector_test::test_3()";

    //
    // Test construction of non-empty matrix with alpha != beta
    //

    ab_orbital_selector s(4, 3);
    if (s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (s.alpha().n_indexes() != 4) {
        fail_test(testname, __FILE__, __LINE__, "# alpha indexes != 4.");
    }
    if (s.beta().n_indexes() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# beta indexes != 3.");
    }

    s.alpha() = orbital_selector(3);
    if (s.alpha().n_indexes() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# alpha indexes != 3.");
    }
}


void ab_orbital_selector_test::test_4a() {

    static const char *testname = "ab_orbital_selector_test::test_4a()";

    //
    // Test change of non-empty matrix from alpha != beta to alpha == beta
    //

    ab_orbital_selector s(4, 6);
    if (s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }

    s.alpha().select(true, 1, 3, 2);
    s.beta().select(false, 2, 5);

    orbital_selector sa = s.alpha();
    s.set_alpha_eq_beta();

    if (! s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta (flag).");
    }
    if (&(s.alpha()) != &(s.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }
    if (s.beta().n_indexes() != 4) {
        fail_test(testname, __FILE__, __LINE__, "# beta elements != 4");
    }
    for (size_t i = 0; i < 4; i++) {
        if (sa.is_selected(true, i) != s.alpha().is_selected(true, i)) {
            fail_test(testname, __FILE__, __LINE__, "alpha changed.");
        }
    }
}


void ab_orbital_selector_test::test_4b() {

    static const char *testname = "ab_orbital_selector_test::test_4b()";

    //
    // Test change of non-empty matrix from alpha == beta to alpha != beta
    //

    size_t n = 5;

    ab_orbital_selector s(n);
    if (! s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha != beta.");
    }

    s.alpha().select(2, 4);
    s.alpha().deselect(3);


    s.set_alpha_neq_beta();
    if (s.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta (flag).");
    }
    if (&(s.alpha()) == &(s.beta())) {
        fail_test(testname, __FILE__, __LINE__, "alpha == beta.");
    }
    if (s.beta().n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# beta selected != 5");
    }
    for (size_t i = 0; i < n; i++) {
        if (s.alpha().is_selected(i) != s.beta().is_selected(i)) {
            fail_test(testname, __FILE__, __LINE__, "alpha != beta (contents)");
        }
    }
}


} // namespace libwfa
