#include <iomanip>
#include <fstream>
#include "ex_ana_printer.h"

namespace libwfa{

using namespace arma;

const double ex_ana_printer::k_au2ang = 0.52917724900;

void ex_ana_printer::perform(ex_analyse &analyse, std::ostream &out) {

    out << std::setw(10) << std::setprecision(6) << std::fixed;
    out << std::endl << std::endl;
    out << std::string(70, '-') << std::endl;
    out << "Performing exciton analysis (TDM)... (in \u212B) " << std::endl;
    out << std::endl;
	out << "For the alpha part:" << std::endl;
	do_print('a', analyse, out);
	if (! analyse.aeqb()) {
	    out << std::endl;
		out << "For the beta part:" << std::endl;
		do_print('b', analyse, out);
	} else {
		out << "Alpha is equal to beta." << std::endl;
	}
	out << std::endl;
	out << std::string(70, '-') << std::endl;
	out << std::endl;

}

void ex_ana_printer::do_print(char spin,
		ex_analyse &analyse, std::ostream &out) {

    out << "<rhx>,<rhy>,<rhz>: "
    		<< analyse.get_rh('x', spin) * k_au2ang << ", "
            << analyse.get_rh('y', spin) * k_au2ang << ", "
            << analyse.get_rh('z', spin) * k_au2ang << std::endl;
    out << "<rex>,<rey>,<rez>: "
    		<< analyse.get_re('x', spin) * k_au2ang << ", "
            << analyse.get_re('y', spin) * k_au2ang << ", "
            << analyse.get_re('z', spin) * k_au2ang << std::endl;
    out << "Electron hole separation |<rh-re>|: "
    		<< analyse.get_sep(spin) * k_au2ang << std::endl;
    out << "Distance for each coord:"
    		<< " x= " << analyse.get_dex_c('x', spin) * k_au2ang
            << " y= " << analyse.get_dex_c('y', spin) * k_au2ang
            << " z= " << analyse.get_dex_c('z', spin) * k_au2ang
            << std::endl;
    out << "Averaged distance over all coordinates: "
            << analyse.get_dex_tot(spin) * k_au2ang << std::endl;
    out << "Sigma (hole): "
    		<< analyse.get_sig_h(spin) * k_au2ang << std::endl;
    out << "Sigma (electron): "
    		<< analyse.get_sig_e(spin) * k_au2ang << std::endl;
    out << "Covariance (rh,re): "
    		<< analyse.get_cov(spin) * k_au2ang * k_au2ang << std::endl;
    out << "Correlation factor (rh,re): "
    		<< analyse.get_corr(spin) << std::endl;
}//end fct


}//end namespace libwfa
