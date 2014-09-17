#include <libwfa/analyses/exciton_analysis_ad.h>
#include <libwfa/core/mom_builder.h>
#include "exciton_analysis_ad_test.h"
#include "exciton_test_data_hhe.h"


namespace libwfa{

void exciton_analysis_ad_test::perform() throw(libtest::test_exception){

    test_1();

}


void exciton_analysis_ad_test::test_1(){

    static const char *testname = "exciton_analysis_ad_test::test_1()";

    try {

    arma::mat s = exciton_test_data_hhe::overlap();
    arma::mat x = exciton_test_data_hhe::x();
    arma::mat y = exciton_test_data_hhe::y();
    arma::mat z = exciton_test_data_hhe::z();
    arma::mat xx = exciton_test_data_hhe::xx();
    arma::mat yy = exciton_test_data_hhe::yy();
    arma::mat zz = exciton_test_data_hhe::zz();

    mom_builder bld(s, x, y, z, xx, yy, zz);

    ab_matrix tdm(true), de(true), at(true);
    tdm.alpha() = exciton_test_data_hhe::tdm();

    // Analyze the particle/hole densities as
    //   attachment/detachment densities
    // Formal negativ sign for the detachment density matrix
    de.alpha() = -tdm.alpha() * s * tdm.alpha().t();
    at.alpha() =  tdm.alpha().t() * s * tdm.alpha();

    exciton_analysis_ad exc(bld, at, de);

    vec rh = exc.moment(false).get(0, 1);
    vec re = exc.moment(false).get(1, 0);
    vec rh2 = exc.moment(false).get(0, 2);
    vec re2 = exc.moment(false).get(2, 0);
    if (norm(rh - exciton_test_data_hhe::rh()) > 1e-4)
        fail_test(testname, __FILE__, __LINE__, "rh");
    if (norm(re - exciton_test_data_hhe::re()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re");
    if (norm(rh2 - exciton_test_data_hhe::rh2()) > 1e-3)
        fail_test(testname, __FILE__, __LINE__, "rh2");
    if (norm(re2 - exciton_test_data_hhe::re2()) > 1e-6)
        fail_test(testname, __FILE__, __LINE__, "re2");

    std::ofstream of("exciton_analysis_ad_test_1", std::ofstream::app);
    exc.analyse(of);

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // end namespace libwfa
