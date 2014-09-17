#include <libwfa/analyses/exciton_analysis.h>
#include <libwfa/core/mom_builder.h>
#include "exciton_test_data_hhe.h"
#include "exciton_analysis_test.h"


namespace libwfa{

void exciton_analysis_test::perform() throw(libtest::test_exception){

    test_1();

}


void exciton_analysis_test::test_1() {

    static const char *testname = "exciton_analysis_test::test_1()";

    try {

    arma::mat s = exciton_test_data_hhe::overlap();
    arma::mat x = exciton_test_data_hhe::x();
    arma::mat y = exciton_test_data_hhe::y();
    arma::mat z = exciton_test_data_hhe::z();
    arma::mat xx = exciton_test_data_hhe::xx();
    arma::mat yy = exciton_test_data_hhe::yy();
    arma::mat zz = exciton_test_data_hhe::zz();

    mom_builder bld(s, x, y, z, xx, yy, zz);

    ab_matrix tdm(true);
    tdm.alpha() = exciton_test_data_hhe::tdm();

    exciton_analysis exc(bld, tdm);

    vec rh = exc.moment(false).get(0, 1);
    vec re = exc.moment(false).get(1, 0);
    vec rh2 = exc.moment(false).get(0, 2);
    vec re2 = exc.moment(false).get(2, 0);
    vec rhre = exc.moment(false).get(1, 1);
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

    std::ofstream of("exciton_analysis_test_1", std::ofstream::app);
    exc.analyse(of);

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // end namespace libwfa
