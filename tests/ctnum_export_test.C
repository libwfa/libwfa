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

    ab_matrix r1(true), r2(true);
    r1.alpha() = symmatu(randu(na, na));
    r2.alpha() = symmatu(randu(na, na));

    ctnum_export("ctnum_export_test_1",
            "Test 1, Set 1 1.23456 0.015").perform(r1);
    ctnum_export("ctnum_export_test_1",
            "Test 1, Set 2, 2.34566, 0.121").perform(r2);
}


void ctnum_export_test::test_2() {

    // Test adjustment of columns

    static const char *testname = "ctnum_export_test::test_2()";

    size_t na = 4;

    ab_matrix r1, r2;
    r1.alpha() = symmatu(randu(na, na));
    r1.beta()  = symmatu(randu(na, na));
    r2.alpha() = symmatu(randu(na + 1, na + 1));
    r2.beta()  = symmatu(randu(na + 1, na + 1));

    ctnum_export("ctnum_export_test_2",
            "Test 2, Set 1, 0.12348, 0.001", 5).perform(r1);
    ctnum_export("ctnum_export_test_2",
            "Test 2, Set 2, 2.34902, 0.382", 5).perform(r2);
}


} // namespace libwfa
