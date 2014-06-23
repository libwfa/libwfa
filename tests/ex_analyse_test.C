#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/core/multipol_con.h>
#include <libwfa/core/transformations_dm.h>
#include "ex_analyse_test.h"
/**TODO: Write Tests for printer. Finish documentation. Look through if layout is okay.
 *
 */

namespace libwfa{

void ex_analyse_test::perform() throw(libtest::test_exception){

    test_ex_total();

} // end fct

using namespace arma;

namespace{

//Test Data Class
class system_data {
private:
    static const double d_tdm[16];
    static const double d_mxx[16];
    static const double d_mx[16];
    static const double d_myy[16];
    static const double d_my[16];
    static const double d_mzz[16];
    static const double d_mz[16];
    static const double d_s[16];

public:
    static Mat<double> overlap() {
        Mat<double> s(d_s,4,4);
        return s;
    }

    static Mat<double> tdm() {
        Mat<double> tdm(d_tdm,4,4);
        return tdm;
    }

    static Mat<double> mxx() {
        Mat<double> mxx(d_mxx,4,4);
        return mxx;
    }

    static Mat<double> mx() {
        Mat<double> mx(d_mx,4,4);
        return mx;
    }

    static Mat<double> myy() {
        Mat<double> myy(d_myy,4,4);
        return myy;
    }

    static Mat<double> my() {
        Mat<double> my(d_my,4,4);
        return my;
    }
    static Mat<double> mzz() {
        Mat<double> mzz(d_mzz,4,4);
        return mzz;
    }

    static Mat<double> mz() {
        Mat<double> mz(d_mz,4,4);
        return mz;
    }

};
//Test Data for H, He case. For further informations, read supplements.
const double system_data::d_tdm[16]={
        0.000000, 0.000000, -0.238998, -0.346058,
        0.000000, 0.000000, -0.459122, -0.664786,
        0.001281, 0.002462, 0.000000, 0.000000,
        -0.001170, -0.002248, 0.000000, 0.000000
};


const double system_data::d_mxx[16]={
        0.261802, 0.307013, 0.000000, 0.000000,
        0.307013, 1.364688, 0.000000, 0.000000,
        0.000000, 0.000000, 0.105856, 0.119238,
        0.000000, 0.000000, 0.119238, 0.652753
};

const double system_data::d_mx[16]={
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000
};

const double system_data::d_myy[16]={
        0.261802, 0.307013, 0.000000, 0.000000,
        0.307013, 1.364688, 0.000000, 0.000000,
        0.000000, 0.000000, 0.105856, 0.119238,
        0.000000, 0.000000, 0.119238, 0.652753
};

const double system_data::d_my[16]={
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000
};

const double system_data::d_mzz[16]={
        0.261802, 0.307013, 0.000000, 0.000000,
        0.307013, 1.364688, 0.000000, 0.000000,
        0.000000, 0.000000, 357.212338, 212.674712,
        0.000000, 0.000000, 212.674712, 357.759236
};

const double system_data::d_mz[16]={
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 18.897261, 11.247951,
        0.000000, 0.000000, 11.247951, 18.897261
};

const double system_data::d_s[16]={
        1.000000, 0.645899, 0.000000, 0.000000,
        0.645899, 1.000000, 0.000000, 0.000000,
        0.000000, 0.000000, 1.000000, 0.595216,
        0.000000, 0.000000, 0.595216, 1.000000
};

}//end unnamed namespaceunit_test

void ex_analyse_test::test_ex_total(){

	static const char *testname = "ex_analyse_test::test_ex_total()";

	arma::Mat<double> s = system_data::overlap();

	ab_matrix tdm(true);
	tdm.alpha()=system_data::tdm();
	arma::Mat<double> mx = system_data::mx();
	arma::Mat<double> my = system_data::my();
	arma::Mat<double> mz = system_data::mz();
	arma::Mat<double> mxx = system_data::mxx();
	arma::Mat<double> myy = system_data::myy();
	arma::Mat<double> mzz = system_data::mzz();

	ab_matrix om(true);
	try{
	form_om(s, tdm, om);
	}catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }


	multipol_con con(mx, mxx, my, myy, mz, mzz, s);

	ex_analyse analyse;
	try{
	analyse.perform(tdm, om, con);
	}catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

	if ((abs(analyse.get_rh('x','a')-0)>1e-6)||
	        (abs(analyse.get_rh('y','a')-0)>1e-6)||
	        (abs(analyse.get_rh('z','a')-18.8972)>1e-6)){
	fail_test(testname, __FILE__, __LINE__,
	        "Exp. Value for h coord. differs to much.");
	}

    if ((abs(analyse.get_rh2('x','a')-0.377513)>1e-6)||
            (abs(analyse.get_rh2('y','a')-0)>1e-6)||
            (abs(analyse.get_rh2('z','a')-357.482)>1e-6)){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for h coord. differs to much.");
    }

    if ((abs(analyse.get_re2('x','a')-0.903178)>1e-6)||
            (abs(analyse.get_re2('y','a')-0)>1e-6)||
            (abs(analyse.get_re2('z','a')-0.904768)>1e-6)){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for h coord. differs to much.");
    }

    if ((abs(analyse.get_rhre('x','a')-0)>1e-6)||
            (abs(analyse.get_rhre('y','a')-0)>1e-6)||
            (abs(analyse.get_rhre('z','a')-0)>1e-6)){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for h coord. differs to much.");
    }

	if ((abs(analyse.get_re('x','a')-0)>1e-6)||
	        (abs(analyse.get_re('y','a')-0)>1e-6)||
	        (abs(analyse.get_re('z','a')-8.41644e-5)>1e-6)){
	fail_test(testname, __FILE__, __LINE__,
	        "Exp. Value for e coord. differs to much.");
	}
	if (abs(analyse.get_sep('a')-18.8971)>1e-6){
	fail_test(testname, __FILE__, __LINE__,
	        "Exp. Value for sep. differs to much.");
	}

    if ((abs(analyse.get_dex_c('x','a')-1.13168)>1e-6)||
            (abs(analyse.get_dex_c('y','a')-0)>1e-6)||
            (abs(analyse.get_dex_c('z','a')-18.9311)>1e-6)){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for distance coord. differs to much.");
    }

    if ((abs(analyse.get_sig_h('a')-0.86985)>1e-6)||
            (abs(analyse.get_sig_e('a')-1.3446)>1e-6)){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for sigma differs to much.");
    }

    if (abs(analyse.get_cov('a')-(-0.00159047))>1e-6){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for cov. differs to much.");
    }

    if (abs(analyse.get_corr('a')-(-0.00135984))>1e-6){
    fail_test(testname, __FILE__, __LINE__,
            "Exp. Value for differs to much. wrong.");
    }



}// end fct

}// end namespace libwfa
