#include "ex_analyse.h"

namespace libwfa {

using namespace arma;


double ex_analyse::ex_multip ( const ab_matrix &tdm, const ab_matrix &op1,
		const ab_matrix &op2, ab_matrix &om, char spin){

	switch (spin){
	case 'a':{
			const Mat<double> &tdm_a=tdm.alpha();
			const Mat<double> &op1_a=op1.alpha();
			const Mat<double> &op2_a=op2.alpha();
			double omval_a=accu(om.alpha());
			double expval_a=0;

			Mat<double> exp_a=(tdm_a*op2_a)%(op1_a*tdm_a);

			expval_a=accu(exp_a)/omval_a;
			return expval_a;
			break;
     }

	case 'b':{
			const Mat<double> &tdm_b=tdm.beta();
			const Mat<double> &op1_b=op1.beta();
			const Mat<double> &op2_b=op2.beta();
			double omval_b=accu(om.beta());
			double expval_b=0;

			Mat<double> exp_b=(tdm_b*op2_b)%(op1_b*tdm_b);

			expval_b=accu(exp_b)/omval_b;

			return expval_b;
			break;
	}
		default:
			return 0;

	} //end switch

}//end fct





}// end namespace libwfa



