#include "ex_analyse_ad.h"

namespace libwfa {

using namespace arma;


void ex_analyse_ad::ex_form_ad(const ab_matrix &att, const ab_matrix &det,
		const contract_i &op){

    m_aeqb = att.is_alpha_eq_beta();
    prom[0] = op.perform(att.alpha(), "s");
    rh[0][0]  = op.perform(det.alpha(), "x")  / prom[0];
    rh[1][0]  = op.perform(det.alpha(), "y")  / prom[0];
    rh[2][0]  = op.perform(det.alpha(), "z")  / prom[0];
    rh2[0][0] = op.perform(det.alpha(), "xx") / prom[0];
    rh2[1][0] = op.perform(det.alpha(), "yy") / prom[0];
    rh2[2][0] = op.perform(det.alpha(), "zz") / prom[0];
    re[0][0]  = op.perform(att.alpha(), "x")  / prom[0];
    re[1][0]  = op.perform(att.alpha(), "y")  / prom[0];
    re[2][0]  = op.perform(att.alpha(), "z")  / prom[0];
    re2[0][0] = op.perform(att.alpha(), "xx") / prom[0];
    re2[1][0] = op.perform(att.alpha(), "yy") / prom[0];
    re2[2][0] = op.perform(att.alpha(), "zz") / prom[0];

    if (!att.is_alpha_eq_beta()){

        prom[1] = op.perform(att.beta(), "s");
        rh[0][1]  = op.perform(det.beta(), "x")  / prom[1];
        rh[1][1]  = op.perform(det.beta(), "y")  / prom[1];
        rh[2][1]  = op.perform(det.beta(), "z")  / prom[1];
        rh2[0][1] = op.perform(det.beta(), "xx") / prom[1];
        rh2[1][1] = op.perform(det.beta(), "yy") / prom[1];
        rh2[2][1] = op.perform(det.beta(), "zz") / prom[1];
        re[0][1]  = op.perform(att.beta(), "x")  / prom[1];
        re[1][1]  = op.perform(att.beta(), "y")  / prom[1];
        re[2][1]  = op.perform(att.beta(), "z")  / prom[1];
        re2[0][1] = op.perform(att.beta(), "xx") / prom[1];
        re2[1][1] = op.perform(att.beta(), "yy") / prom[1];
        re2[2][1] = op.perform(att.beta(), "zz") / prom[1];
    } else {
    	prom[1] = prom[0];
        rh[0][1]  = rh[0][0];
        rh[1][1]  = rh[1][0];
        rh[2][1]  = rh[2][0];
        rh2[0][1] = rh2[0][0];
        rh2[1][1] = rh2[1][0];
        rh2[2][1] = rh2[2][0];
        re[0][1]  = re[0][0];
        re[1][1]  = re[1][0];
        re[2][1]  = re[2][0];
        re2[0][1] = re2[0][0];
        re2[1][1] = re2[1][0];
        re2[2][1] = re2[2][0];
    }//end else
}//end fct

double ex_analyse_ad::ex_mean_sep_ad (char spin){

	size_t s = determine_spin(spin);

    double sum = 0;
    for (size_t i = 0; i < 3; i++) {
    	double tmp = rh[i][s] - re[i][s];
    	sum += tmp * tmp;
    }
    return sqrt(sum);
} //end fct

double ex_analyse_ad::ex_sig_h_ad (char spin){

	size_t s = determine_spin(spin);
	double sum = 0;
	for (size_t i = 0; i < 3; i++) {
		sum += rh2[i][s] - rh[i][s] * rh[i][s];
	}
	return sqrt(sum);
}//end fct

double ex_analyse_ad::ex_sig_e_ad (char spin){

	size_t s = determine_spin(spin);
	double sum = 0;
	for (size_t i = 0; i < 3; i++) {
		sum += re2[i][s] - re[i][s] * re[i][s];
	}
	return sqrt(sum);
}//end fct

void ex_analyse_ad::perform(const ab_matrix &att, const ab_matrix &det,
	const contract_i &op){
//Forming all needed values for the calculations;
    ex_form_ad(att, det, op);
//Performing all calculations, saving them in the resp. arrays.

    sep[0]=ex_mean_sep_ad('a');
    sig_h[0]=ex_sig_h_ad('a');
    sig_e[0]=ex_sig_e_ad('a');


    if (!att.is_alpha_eq_beta()){ //If beta is not equal alpha

        sep[1] = ex_mean_sep_ad('b');
        sig_h[1] = ex_sig_h_ad('b');
        sig_e[1] = ex_sig_e_ad('b');

    }else{

        sep[1] = sep[0];
        sig_h[1] = sig_h[0];
        sig_e[1] = sig_e[0];

    }//end else
}//end fct

double ex_analyse_ad::get_sig_h(char spin) {

	return sig_h[determine_spin(spin)];
} //end fct

double ex_analyse_ad::get_sig_e(char spin) {

	return sig_e[determine_spin(spin)];
}//end fct

double ex_analyse_ad::get_sep(char spin) {

	return sep[determine_spin(spin)];
} //end fct

double ex_analyse_ad::get_rh (char koord, char spin){

	return rh[determine_coord(koord)][determine_spin(spin)];
}//end fct

double ex_analyse_ad::get_re (char koord, char spin){

	return re[determine_coord(koord)][determine_spin(spin)];
}//end fct

double ex_analyse_ad::get_rh2 (char koord, char spin){

	return rh2[determine_coord(koord)][determine_spin(spin)];
}//end fct

double ex_analyse_ad::get_re2 (char koord, char spin) {

	return re2[determine_coord(koord)][determine_spin(spin)];
}//end fct


size_t ex_analyse_ad::determine_coord(char coord) {

	if (coord == 'x') return 0;
	else if (coord == 'y') return 1;
	else if (coord == 'z') return 2;
	else
		throw 1;
}//endfct


size_t ex_analyse_ad::determine_spin(char spin) {

	if (spin == 'a') return 0;
	else if (spin == 'b') return 1;
	else
		throw 1;
}//endfct

} //end namespace
