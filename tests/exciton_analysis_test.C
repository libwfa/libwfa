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

    ab_exciton_moments mom;

    try {

    arma::Mat<double> s = exciton_test_data_hhe::overlap();
    arma::Mat<double> x = exciton_test_data_hhe::x();
    arma::Mat<double> y = exciton_test_data_hhe::y();
    arma::Mat<double> z = exciton_test_data_hhe::z();
    arma::Mat<double> xx = exciton_test_data_hhe::xx();
    arma::Mat<double> yy = exciton_test_data_hhe::yy();
    arma::Mat<double> zz = exciton_test_data_hhe::zz();

    mom_builder bld(s, x, y, z, xx, yy, zz);

    ab_matrix tdm(true);
    tdm.alpha() = exciton_test_data_hhe::tdm();

    exciton_analysis(bld, tdm).perform(mom);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    Col<double> rh = mom.alpha().get(0, 1);
    if (norm(rh - exciton_test_data_hhe::rh()) > 1e-4) {
        fail_test(testname, __FILE__, __LINE__, "rh");
    }
    Col<double> re = mom.alpha().get(1, 0);
    if (norm(re - exciton_test_data_hhe::re()) > 1e-6) {
        fail_test(testname, __FILE__, __LINE__, "re");
    }
    Col<double> rh2 = mom.alpha().get(0, 2);
    if (norm(rh2 - exciton_test_data_hhe::rh2()) > 1e-3) {
        fail_test(testname, __FILE__, __LINE__, "rh2");
    }
    Col<double> re2 = mom.alpha().get(2, 0);
    if (norm(re2 - exciton_test_data_hhe::re2()) > 1e-6) {
        fail_test(testname, __FILE__, __LINE__, "re2");
    }
    Col<double> rhre = mom.alpha().get(1, 1);
    if (norm(rhre - exciton_test_data_hhe::rhre()) > 1e-6) {
        fail_test(testname, __FILE__, __LINE__, "rhre");
    }
}


} // end namespace libwfa
