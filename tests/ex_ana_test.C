#include <libwfa/core/transformations_dm.h>
#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/export/ex_ana_printer.h>
#include <libwfa/analyses/ex_analyse_ad.h>
#include <libwfa/export/ex_ana_printer_ad.h>
#include "ex_ana_test.h"


namespace libwfa{

void ex_ana_test::perform() throw(libtest::test_exception){

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

void ex_ana_test::test_ex_total(){
    static const char *testname = "ex_ana_test::test_ex_total()";
    //Forming ab_matrix for test purposes, where alpha==beta
        ab_matrix tdm(4, 4);
        tdm.alpha()=system_data::tdm();
        tdm.set_alpha_eq_beta();

        ab_matrix mxx(4, 4);
        mxx.alpha()=system_data::mxx();
        mxx.set_alpha_eq_beta();

        ab_matrix mx(4, 4);
        mx.alpha()=system_data::mx();
        mx.set_alpha_eq_beta();

        ab_matrix myy(4, 4);
        myy.alpha()=system_data::myy();
        myy.set_alpha_eq_beta();

        ab_matrix my(4, 4);
        my.alpha()=system_data::my();
        my.set_alpha_eq_beta();

        ab_matrix mzz(4, 4);
        mzz.alpha()=system_data::mzz();
        mzz.set_alpha_eq_beta();

        ab_matrix mz(4, 4);
        mz.alpha()=system_data::mz();
        mz.set_alpha_eq_beta();

        ab_matrix s(4, 4);
        s.alpha()=system_data::overlap();
        s.set_alpha_eq_beta();

        ab_matrix om(4, 4);
        form_om (system_data::overlap(), tdm, om);
        om.set_alpha_eq_beta();

        //Creating the test and printer object, performing all tests and
        //printing the results
        ex_analyse analyse;
        ex_ana_printer ana_p;
        analyse.perform(tdm,s,mxx,mx,myy,my,mzz,mz,om);
        ana_p.perform(tdm.is_alpha_eq_beta(), analyse, cout);

        ex_analyse_ad analyse_ad;
        ex_ana_printer_ad ana_p_ad;

        //Creating Detachment, Attachment Matrices
        ab_matrix det(4,4);
        det.alpha()=tdm.alpha()*s.alpha()*tdm.alpha().t();
        det.set_alpha_eq_beta();

        ab_matrix att(4,4);
        att.alpha()=tdm.alpha().t()*s.alpha()*tdm.alpha();
        att.set_alpha_eq_beta();

        analyse_ad.perform(att,det,mx,mxx,my,myy,mz,mzz,om);
        ana_p_ad.perform(tdm.is_alpha_eq_beta(), analyse_ad, cout);



}// end fct

}// end namespace libwfa
