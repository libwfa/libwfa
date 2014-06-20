#include "ex_analyse.h"

namespace libwfa {

using namespace arma;



void ex_analyse::ex_form (const ab_matrix &tdm, const ab_matrix &om,
        const contract_i &op){

	m_aeqb = tdm.is_alpha_eq_beta();
	double om_tot = accu(om.alpha());

    rh[0][0] = op.perform(tdm.alpha(), "s", "x") / om_tot;
    rh[1][0] = op.perform(tdm.alpha(), "s", "y") / om_tot;
    rh[2][0] = op.perform(tdm.alpha(), "s", "z") / om_tot;

    re[0][0] = op.perform(tdm.alpha(), "x", "s") / om_tot;
    re[1][0] = op.perform(tdm.alpha(), "y", "s") / om_tot;
    re[2][0] = op.perform(tdm.alpha(), "z", "s") / om_tot;

    rh2[0][0] = op.perform(tdm.alpha(), "s", "xx") / om_tot;
    rh2[1][0] = op.perform(tdm.alpha(), "s", "yy") / om_tot;
    rh2[2][0] = op.perform(tdm.alpha(), "s", "zz") / om_tot;

    re2[0][0] = op.perform(tdm.alpha(), "xx", "s") / om_tot;
    re2[1][0] = op.perform(tdm.alpha(), "yy", "s") / om_tot;
    re2[2][0] = op.perform(tdm.alpha(), "zz", "s") / om_tot;

    rhre[0][0] = op.perform(tdm.alpha(), "x", "x") / om_tot;
    rhre[1][0] = op.perform(tdm.alpha(), "y", "y") / om_tot;
    rhre[2][0] = op.perform(tdm.alpha(), "z", "z") / om_tot;

    if (!tdm.is_alpha_eq_beta()){

    	om_tot = accu(om.beta());

        rh[0][1] = op.perform(tdm.beta(), "s", "x") / om_tot;
        rh[1][1] = op.perform(tdm.beta(), "s", "y") / om_tot;
        rh[2][1] = op.perform(tdm.beta(), "s", "z") / om_tot;

        re[0][1] = op.perform(tdm.beta(), "x", "s") / om_tot;
        re[1][1] = op.perform(tdm.beta(), "y", "s") / om_tot;
        re[2][1] = op.perform(tdm.beta(), "z", "s") / om_tot;

        rh2[0][1] = op.perform(tdm.beta(), "s", "xx") / om_tot;
        rh2[1][1] = op.perform(tdm.beta(), "s", "yy") / om_tot;
        rh2[2][1] = op.perform(tdm.beta(), "s", "zz") / om_tot;

        re2[0][1] = op.perform(tdm.beta(), "xx", "s") / om_tot;
        re2[1][1] = op.perform(tdm.beta(), "yy", "s") / om_tot;
        re2[2][1] = op.perform(tdm.beta(), "zz", "s") / om_tot;

        rhre[0][1] = op.perform(tdm.beta(), "x", "x") / om_tot;
        rhre[1][1] = op.perform(tdm.beta(), "y", "y") / om_tot;
        rhre[2][1] = op.perform(tdm.beta(), "z", "z") / om_tot;
    }
    else{
        rh[0][1]=rh[0][0];
        rh[1][1]=rh[1][0];
        rh[2][1]=rh[2][0];
        rh2[0][1]=rh2[0][0];
        rh2[1][1]=rh2[1][0];
        rh2[2][1]=rh2[2][0];
        re[0][1]=re[0][0];
        re[1][1]=re[1][0];
        re[2][1]=re[2][0];
        re2[0][1]=re2[0][0];
        re2[1][1]=re2[1][0];
        re2[2][1]=re2[2][0];
        rhre[0][1]=rhre[0][0];
        rhre[1][1]=rhre[1][0];
        rhre[2][1]=rhre[2][0];
    }//end else

} // end fct

double ex_analyse::ex_d_ex_c (char koord, char spin){

    size_t c = determine_coord(koord);
    size_t s = determine_spin(spin);

    return sqrt(rh2[c][s] - 2 * rhre[c][s] + re2[c][s]);
}//end fct

double ex_analyse::ex_d_ex_tot (char spin){
//TODO:discuss which would be better, private=will only be called by the class, additional calc not needed
	return sqrt(ex_d_ex_c('x', spin) * ex_d_ex_c('x', spin)
			+ ex_d_ex_c('y', spin) * ex_d_ex_c('y', spin)
			+ ex_d_ex_c('z', spin) * ex_d_ex_c('z', spin));
	/**
	 size_t s=determine_spin(spin);
	 return sqrt(dex_c[determine_coord('x')][s]
	             * dex_c[determine_coord('x')][s]
	             + dex_c[determine_coord('y')][s]
	             * dex_c[determine_coord('y')][s]
	             + dex_c[determine_coord('z')][s]
	             * dex_c[determine_coord('z')][s]);
	 **/

}//end fct

double ex_analyse::get_rh (char koord, char spin){

	return rh[determine_coord(koord)][determine_spin(spin)];
}

double ex_analyse::get_re (char koord, char spin){

	return re[determine_coord(koord)][determine_spin(spin)];
}

double ex_analyse::get_rh2 (char koord, char spin){

	return rh2[determine_coord(koord)][determine_spin(spin)];
}

double ex_analyse::get_re2 (char koord, char spin){

	return re2[determine_coord(koord)][determine_spin(spin)];
}

double ex_analyse::get_rhre (char koord, char spin){

	return rhre[determine_coord(koord)][determine_spin(spin)];
}

double ex_analyse::ex_mean_sep (char spin){

	size_t s = determine_spin(spin);
	double sum = 0;
	for (size_t i = 0; i < 3; i++) {
		double tmp = rh[i][s] - re[i][s];
		sum += tmp * tmp;
	}
	return sqrt(sum);
} //end fct

double ex_analyse::ex_sig_h (char spin){

	size_t s = determine_spin(spin);
	double sum = 0;
	for (size_t i = 0; i < 3; i++) sum += (rh2[i][s] - rh[i][s] * rh[i][s]);
	return sqrt(sum);
}//end fct

double ex_analyse::ex_sig_e (char spin){

	size_t s = determine_spin(spin);
	double sum = 0;
	for (size_t i = 0; i < 3; i++) sum += (re2[i][s] - re[i][s] * re[i][s]);
	return sqrt(sum);
}//end fct

double ex_analyse::ex_cov (char spin){

	size_t s = determine_spin(spin);
	double erg = 0;
	for (size_t i = 0; i < 3; i++) erg += (rhre[i][s] - rh[i][s] * re[i][s]);
	return erg;

}//end fct

double ex_analyse::ex_corr(char spin){

//TODO:discuss which better, private=will only be called by the class, additional calc not needed

    return ex_cov(spin) / (ex_sig_e(spin) * ex_sig_h(spin));
    /**
     size_t s=determine_spin(spin);
     return cov[s] / (sig_e[s]*sig_h[s]);
     */
}


double ex_analyse::get_sep(char spin) {

	return sep[determine_spin(spin)];
} //end fct

double ex_analyse::get_dex_c(char coord, char spin) {

	return dex_c[determine_coord(coord)][determine_spin(spin)];
} //end fct

double ex_analyse::get_dex_tot(char spin) {

	return dex_tot[determine_spin(spin)];
} //end fct

double ex_analyse::get_sig_h(char spin) {

	return sig_h[determine_spin(spin)];
} //end fct

double ex_analyse::get_sig_e(char spin) {

	return sig_e[determine_spin(spin)];

}//end fct

double ex_analyse::get_cov(char spin) {

	return cov[determine_spin(spin)];
} //end fct

double ex_analyse::get_corr(char spin) {

	return corr[determine_spin(spin)];
} //end fct


void ex_analyse::perform(const ab_matrix &tdm, const ab_matrix &om,
        const contract_i &name){
//Forming all needed values for the calculations;
    ex_form(tdm, om, name);
//Performing all calculations, saving them in the resp. arrays.

    sep[0]=ex_mean_sep('a');
    dex_c[0][0]=ex_d_ex_c ('x', 'a');
    dex_c[1][0]=ex_d_ex_c ('y', 'a');
    dex_c[2][0]=ex_d_ex_c ('z', 'a');
    dex_tot[0]=ex_d_ex_tot('a');
    sig_h[0]=ex_sig_h('a');
    sig_e[0]=ex_sig_e('a');
    cov[0]=ex_cov('a');
    corr[0]=ex_corr('a');

    if (!tdm.is_alpha_eq_beta()){ //If beta is not equal alpha

        sep[1] = ex_mean_sep('b');
        dex_c[0][1] = ex_d_ex_c('x', 'b');
        dex_c[1][1] = ex_d_ex_c('y', 'b');
        dex_c[2][1] = ex_d_ex_c('z', 'b');
        dex_tot[1] = ex_d_ex_tot('b');
        sig_h[1] = ex_sig_h('b');
        sig_e[1] = ex_sig_e('b');
        cov[1] = ex_cov('b');
        corr[1] = ex_corr('b');

    }else{

        sep[1] = sep[0];
        dex_c[0][1] = dex_c[0][0];
        dex_c[1][1] = dex_c[1][0];
        dex_c[2][1] = dex_c[2][0];
        dex_tot[1] = dex_tot[0];
        sig_h[1] = sig_h[0];
        sig_e[1] = sig_e[0];
        cov[1] = cov[0];
        corr[1] = corr[0];

    }//end else
}//end fct


size_t ex_analyse::determine_coord(char coord) {

	if (coord == 'x') return 0;
	else if (coord == 'y') return 1;
	else if (coord == 'z') return 2;
	else
		throw 1;
}


size_t ex_analyse::determine_spin(char spin) {

	if (spin == 'a') return 0;
	else if (spin == 'b') return 1;
	else
		throw 1;
}

}// end namespace libwfa





