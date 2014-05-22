#include <libwfa/core/transformations_dm.h>
#include <libwfa/analyses/ex_analyse.h>
#include "ex_ana_test.h"

namespace libwfa{

void ex_ana_test::perform() throw(libtest::test_exception){

	test_ex_multip();

} // end fct

using namespace arma;

namespace{

//Test Data
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

const double system_data::d_tdm[16]={
		0.000000, 0.000000,	-0.238998, -0.346058,
		0.000000, 0.000000, -0.459122, -0.664786,
		0.001281, 0.002462, 0.000000, 0.000000,
		-0.001170, -0.002248, 0.000000,	0.000000
};


const double system_data::d_mxx[16]={
		0.261802, 0.307013,	0.000000, 0.000000,
		0.307013, 1.364688,	0.000000, 0.000000,
		0.000000, 0.000000,	0.105856, 0.119238,
		0.000000, 0.000000,	0.119238, 0.652753
};

const double system_data::d_mx[16]={
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000
};

const double system_data::d_myy[16]={
		0.261802, 0.307013,	0.000000, 0.000000,
		0.307013, 1.364688,	0.000000, 0.000000,
		0.000000, 0.000000,	0.105856, 0.119238,
		0.000000, 0.000000,	0.119238, 0.652753
};

const double system_data::d_my[16]={
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000
};

const double system_data::d_mzz[16]={
		0.261802, 0.307013,	0.000000, 0.000000,
		0.307013, 1.364688,	0.000000, 0.000000,
		0.000000, 0.000000,	357.212338, 212.674712,
		0.000000, 0.000000,	212.674712, 357.759236
};

const double system_data::d_mz[16]={
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	18.897261, 11.247951,
		0.000000, 0.000000,	11.247951, 18.897261
};

const double system_data::d_s[16]={
		1.000000, 0.645899,	0.000000, 0.000000,
		0.645899, 1.000000,	0.000000, 0.000000,
		0.000000, 0.000000,	1.000000, 0.595216,
		0.000000, 0.000000,	0.595216, 1.000000
};

}//end unnamed namespace

void ex_ana_test::test_ex_multip(){
	static const char *testname = "ex_ana_test::test_ex_multip()";

	//
	// Test calculating the expected value of the averaged distance for an exciton over all spacial coordinates.
	//
	ex_analyse analyse;

	// Forming the reference matrices out of the given values, setting alpha==beta.
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

	//Setting the given value for the pre-calculated result. Setting all newly generated var to zero.
	double erg_x=0;
	double erg_y=0;
	double erg_z=0;
	double erg=0;
	double erg_ref=10.05365128; //Pre-Calculated in excel

	//Trying to calculate the results for d_x²,d_y²,d_z²; if it doesn't work, catching the error and giving out a message.
	try{
		erg_x=analyse.ex_multip(tdm, mxx,s,om,'a')-2*analyse.ex_multip(tdm,mx,mx,om,'a')+analyse.ex_multip(tdm,s,mxx,om,'a');
		}catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
		}

	try{
		erg_y=analyse.ex_multip(tdm, myy,s,om,'a')-2*analyse.ex_multip(tdm,my,my,om,'a')+analyse.ex_multip(tdm,s,myy,om,'a');
		}catch(std::exception &e) {
		fail_test(testname, __FILE__, __LINE__, e.what());
		}

	try{
		erg_z=analyse.ex_multip(tdm, mzz,s,om,'a')-2*analyse.ex_multip(tdm,mz,mz,om,'a')+analyse.ex_multip(tdm,s,mzz,om,'a');
		}catch(std::exception &e) {
		fail_test(testname, __FILE__, __LINE__, e.what());
		}

	erg=sqrt(erg_x+erg_y+erg_z)*0.52917724900; //calculating the final d_av over all spacial coordinates, multiplying with a conversion value.
	//Testing if the calculated value differs too much from the pre-calculated one.
	//std::cout<<"Diff: "<<std::abs(erg-erg_ref)<<" Omega: "<<accu(om.alpha())<<" X:"<<erg_x<<" Y: "<<erg_y<<" Z: "<<erg_z<<" The result is: "<<erg<<" The reference is: "<<erg_ref<<" ";
	std::cout<<"The difference between the calculated and the reference value is: "<<std::abs(erg-erg_ref)<<" ";
	if (std::abs(erg-erg_ref)>1e-6){
	    fail_test(testname, __FILE__, __LINE__, "The results do not match.");
	}
	else  {
		std::cout<< "The results do match. ";
	}

}// end fct

}// end namespace libwfa
