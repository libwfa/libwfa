#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/ctnum_export.h>
#include "ctnum_export_test.h"

namespace libwfa {

using namespace arma;


void ctnum_export_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
}


void ctnum_export_test::test_1() {

    // Basic functionality test

    static const char *testname = "ctnum_export_test::test_1()";

    size_t na = 4;

    double tot[2];
    ab_matrix r1(true), r2(true);
    r1.alpha() = symmatu(randu(na, na));
    r2.alpha() = symmatu(randu(na, na));

    std::ostringstream ss1;
    ctnum_export pr1, pr2;
    pr1.set_state_info("test_1_set_1", "Test 1, Set 1", 1.23456, 0.015);
    pr1.perform(r1, tot, ss1);
    pr2.set_state_info("test_1_set_2", "Test 1, Set 2", 2.34566, 0.121);
    pr2.perform(r2, tot, ss1);
}


void ctnum_export_test::test_2() {

    // Test adjustment of columns

    static const char *testname = "ctnum_export_test::test_2()";

    size_t na = 4;

    double tot[2];
    ab_matrix r1, r2;
    r1.alpha() = symmatu(randu(na, na));
    r1.beta()  = symmatu(randu(na, na));
    r2.alpha() = symmatu(randu(na + 1, na + 1));
    r2.beta()  = symmatu(randu(na + 1, na + 1));

    std::ostringstream ss2;
    ctnum_export pr1(5), pr2(5);
    pr1.set_state_info("test_2_set_1", "Test 2, Set 1", 0.12348, 0.001);
    pr1.perform(r1, tot, ss2);
    pr2.set_state_info("test_2_set_2", "Test 2, Set 2", 2.34902, 0.382);
    pr2.perform(r2, tot, ss2);
}


} // namespace libwfa
