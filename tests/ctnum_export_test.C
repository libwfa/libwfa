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

    Mat<double> r(na, na);
    r.randu();

    ctnum_data ct;
    ct.add("set_1", 1.23456, 0.015) = symmatu(r);
    ct.add("set_2", 2.34566, 0.121) = symmatl(r);

    ctnum_export pr("test1_");
    pr.perform(ct);

}


void ctnum_export_test::test_2() {

    // Test adjustment of columns

    static const char *testname = "ctnum_export_test::test_2()";

    size_t na = 4;

    Mat<double> r1(na, na), r2(na + 1, na + 1);
    r1.randu(); r2.randu();

    ctnum_data ct;
    ct.add("set_1", 0.12348, 0.001) = symmatu(r1);
    ct.add("set_2", 1.30566, 0.012) = symmatl(r1);
    ct.add("set_3", 2.34902, 0.382) = symmatu(r2);
    ct.add("set_4", 2.93426, 0.102) = symmatl(r2);

    ctnum_export pr("test2_", 5);
    pr.perform(ct);
}


} // namespace libwfa
