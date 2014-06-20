#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/core/contract.h>
#include <libwfa/core/transformations_dm.h>
#include "ex_analyse_test.h"


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

	//Forming ab_matrix for test purposes, where alpha==beta
	arma::Mat<double> s = system_data::overlap();

	ab_matrix tdm(true);
	arma::Mat<double> mx = system_data::mx();
	arma::Mat<double> my = system_data::my();
	arma::Mat<double> mz = system_data::mz();
	arma::Mat<double> mxx = system_data::mxx();
	arma::Mat<double> myy = system_data::myy();
	arma::Mat<double> mzz = system_data::mzz();

	ab_matrix om(true);
	form_om(s, tdm, om);

	contract con(mx, mxx, my, myy, mz, mzz, s);

	ex_analyse analyse;

	analyse.perform(tdm, om, con);

	// TODO: compare results in ex_analyse to reference data

	// The idea of a unit test is to provide an easy control on the correctness
	// of a class / function. Whenever someone changes the code inside of
	// a function / class, running the unit test should report, if the new
	// code is working OK. Thus, writing a unit test which does not do an
	// error check in the end is almost completely useless. There are only
	// a few exceptions to this:
	// 1) The test should check on the error handling of a function / class.
	//    In this case one would probably only catch the exception thrown by
	//    the function / class and either forward it or accept the test as
	//    working (depending on what's being tested).
	// 2) The function / class is exclusively for printing. Even then, if
	//    the world was perfect and the output stream were available as string,
	//    the result should be analysed. However, analysing a string is far more
	//    effort than just looking at it.

	// In this case however (with the printer being separated from the analysis)
	// you can easily compare the results of the analysis numerically to
	// reference data. As example how you could do this look at
	// void transformations_dm_test::test_form_eh_1b().

	fail_test(testname, __FILE__, __LINE__, "No check");

}// end fct

}// end namespace libwfa
