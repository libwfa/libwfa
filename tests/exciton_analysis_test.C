#include <libwfa/analyses/exciton_analysis.h>
#include <libwfa/core/mom_builder.h>
#include "exciton_test_data_hhe.h"
#include "exciton_analysis_test.h"


namespace libwfa{

void exciton_analysis_test::perform() throw(libtest::test_exception){

    test_1();
    test_2();
}


using namespace arma;

// The tests have the same data.
// Only difference: test_2 is unrestricted.

void exciton_analysis_test::test_1() {

    static const char *testname = "exciton_analysis_test::test_1()";

    try {

    mat s = exciton_test_data_hhe::overlap();
    mat x = exciton_test_data_hhe::x();
    mat y = exciton_test_data_hhe::y();
    mat z = exciton_test_data_hhe::z();
    mat xx = exciton_test_data_hhe::xx();
    mat yy = exciton_test_data_hhe::yy();
    mat zz = exciton_test_data_hhe::zz();

    mom_builder bld(s, x, y, z, xx, yy, zz);

    ab_matrix tdm(true);
    tdm.alpha() = exciton_test_data_hhe::tdm();

    exciton_analysis exc(bld, tdm);

    vec rh = exc.moment(false).get(0, 1);
    vec re = exc.moment(false).get(1, 0);
    vec rh2 = exc.moment(false).get(0, 2);
    vec re2 = exc.moment(false).get(2, 0);
    vec rhre = exc.moment(false).get(1, 1);
    vec tdip = exc.moment(false).get(1);
    vec tquad = exc.moment(true).get(2);
    if (norm(rh - exciton_test_data_hhe::rh()) > 1e-4)
        fail_test(testname, __FILE__, __LINE__, "rh");
    if (norm(re - exciton_test_data_hhe::re()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re");
    if (norm(rh2 - exciton_test_data_hhe::rh2()) > 1e-3)
        fail_test(testname, __FILE__, __LINE__, "rh2");
    if (norm(re2 - exciton_test_data_hhe::re2()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re2");
    if (norm(rhre - exciton_test_data_hhe::rhre()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "rhre");
    if (norm(tdip - exciton_test_data_hhe::tdip()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "tdip");
    if (norm(tquad - exciton_test_data_hhe::tquad()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "tquad");

    std::ofstream of("exciton_analysis_test_1.out", std::ofstream::app);
    exc.analyse(of);

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}

void exciton_analysis_test::test_2() {

    static const char *testname = "exciton_analysis_test::test_2()";

    try {

    mat s = exciton_test_data_hhe::overlap();
    mat x = exciton_test_data_hhe::x();
    mat y = exciton_test_data_hhe::y();
    mat z = exciton_test_data_hhe::z();
    mat xx = exciton_test_data_hhe::xx();
    mat yy = exciton_test_data_hhe::yy();
    mat zz = exciton_test_data_hhe::zz();

    mom_builder bld(s, x, y, z, xx, yy, zz);

    ab_matrix tdm(false);
    tdm.alpha() = exciton_test_data_hhe::tdm();
    tdm.beta()  = exciton_test_data_hhe::tdm();

    exciton_analysis exc(bld, tdm);

    vec rh = exc.moment(true).get(0, 1);
    vec re = exc.moment(true).get(1, 0);
    vec rh2 = exc.moment(true).get(0, 2);
    vec re2 = exc.moment(true).get(2, 0);
    vec rhre = exc.moment(true).get(1, 1);
    vec tdip = exc.moment(true).get(1);
    vec tquad = exc.moment(true).get(2);
    if (norm(rh - exciton_test_data_hhe::rh()) > 1e-4)
        fail_test(testname, __FILE__, __LINE__, "rh");
    if (norm(re - exciton_test_data_hhe::re()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re");
    if (norm(rh2 - exciton_test_data_hhe::rh2()) > 1e-3)
        fail_test(testname, __FILE__, __LINE__, "rh2");
    if (norm(re2 - exciton_test_data_hhe::re2()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re2");
// Tests for transition dipole and qquadrupole
// -> Divide by two to get alpha/beta component
    if (norm(tdip - 0.5*exciton_test_data_hhe::tdip()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "tdip");
    if (norm(tquad - 0.5*exciton_test_data_hhe::tquad()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "tquad");


    std::ofstream of("exciton_analysis_test_2.out", std::ofstream::app);
    exc.analyse(of);

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}

} // end namespace libwfa
