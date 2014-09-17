#include <libwfa/core/mom_builder.h>
#include "mom_builder_test.h"

namespace libwfa {

using namespace arma;


void mom_builder_test::perform() throw(libtest::test_exception) {

    test_1();
}


void mom_builder_test::test_1() {

    static const char *testname = "mom_builder_test::test_1()";

    try {

    size_t nao = 4;

    mat s =  randu< mat >(nao, nao); s = s + s.t();
    mat x =  randu< mat >(nao, nao); x = x + x.t();
    mat y =  randu< mat >(nao, nao); y = y + y.t();
    mat z =  randu< mat >(nao, nao); z = z + z.t();
    mat xx = randu< mat >(nao, nao); xx = xx + xx.t();
    mat yy = randu< mat >(nao, nao); yy = yy + yy.t();
    mat zz = randu< mat >(nao, nao); zz = zz + zz.t();
    mat dm = randu< mat >(nao, nao); dm = dm + dm.t();

    mom_builder m(s, x, y, z, xx, yy, zz);
    double a = m.perform(dm, 0, 1, 1, 1);
    double b = m.perform(dm, 0, 0);

    double a_ref = trace(dm * x * dm * y);
    double b_ref = trace(dm * s);

    if (fabs(a - a_ref) > 1e-12) {
        fail_test(testname, __FILE__, __LINE__, "a");
    }
    if (fabs(b - b_ref) > 1e-12) {
        fail_test(testname, __FILE__, __LINE__, "b");
    }

    } catch (libtest::test_exception &e) {
        throw;
    } catch (std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
