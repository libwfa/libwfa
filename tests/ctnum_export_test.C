#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/ctnum_export.h>
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
    ctnum_export pr1(ss1, "test1_set_1"), pr2(ss1, "test_1_set_2");
    pr1.set_state_info(1.23456, 0.015);
    pr1.perform(r1, tot);
    pr2.set_state_info(2.34566, 0.121);
    pr2.perform(r2, tot);


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
    ctnum_export pr1(ss2, "test_2_set_1", 5), pr2(ss2, "test_2_set_2", 5);
    pr1.set_state_info(0.12348, 0.001);
    pr1.perform(r1, tot);
    pr2.set_state_info(2.34902, 0.382);
    pr2.perform(r2, tot);
}


} // namespace libwfa
