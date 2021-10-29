#include <libwfa/analyses/exciton_analysis.h>
#include <libwfa/core/mom_builder.h>
#include "exciton_test_data_hhe.h"
#include "non_orth_test.h"


namespace libwfa{

void non_orth_test::perform() throw(libtest::test_exception){

    test_1();
    //test_2();
}


using namespace arma;

void non_orth_test::test_1() {

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
    std::ofstream of("non_orth_test_1", std::ofstream::app);

    // Create dummy MO-coefficients
    ab_matrix c_ref(false), c_exc(false);
    {
        mat u;
        vec e;
        arma::eig_sym(e, u, s);    
        
        // Get reference MO-coefficients via a Lowdin orthogonalisation
        c_ref.alpha() = u * diagmat(1./sqrt(e)) * u.t();
        c_ref.beta() = u * diagmat(1./sqrt(e)) * u.t();

        c_ref.alpha().print(of, "c_ref");
        of << "This should be a unit matrix" << std::endl;
        (c_ref.alpha() * s * c_ref.alpha().t()).print(of, "CSC^T");

        // TODO: construct the "MOM" MO-coefficients via a unitary transformation of c_ref
    }
    
    ab_matrix tdm(false);
    // TODO: construct the 1TDM via Andrew Gilbert's formula
    tdm.alpha() = exciton_test_data_hhe::tdm();
    tdm.beta()  = exciton_test_data_hhe::tdm();

    exciton_analysis exc(bld, tdm);

    exc.analyse(of);

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}

} // end namespace libwfa
