#include <iomanip>
#include <fstream>
#include "ex_ana_printer.h"

namespace libwfa{

using namespace arma;

void ex_ana_printer::perform(bool aeqb, ex_analyse &analyse, std::ostream &out){

    do_print (aeqb, analyse, out);

}//end fct

void ex_ana_printer::do_print(bool aeqb, ex_analyse &analyse, std::ostream &out){
    double c=0.52917724900;
    out<< std::setw(10) << std::setprecision(6) << std::fixed;
    out<<"\n \n";
    out<<"---------------------------------------------------------------------"
            << endl;
    out << "Performing exciton analysis (TDM)... (in \u212B) " << endl;

    out << "\n";
    out << "For the alpha part:" << endl;
    out << "<rhx>,<rhy>,<rhz>: " << analyse.get_rh('x', 'a') * c << ", "
            << analyse.get_rh('y', 'a') * c << ", "
            << analyse.get_rh('z', 'a') * c << endl;
    out << "<rex>,<rey>,<rez>: " << analyse.get_re('x', 'a') * c << ", "
            << analyse.get_re('y', 'a') * c << ", "
            << analyse.get_re('z', 'a') * c << endl;
    out << "Electron hole separation |<rh-re>|: " << analyse.get_sep('a') * c
            << endl;
    out << "Distance for each coord: x= " << analyse.get_dex_c('x', 'a') * c
            << " y= " << analyse.get_dex_c('y', 'a') * c << " z="
            << analyse.get_dex_c('z', 'a') * c << endl;
    out << "Averaged distance over all coordinates: "
            << analyse.get_dex_tot('a') * c << endl;
    out << "Sigma (hole): " << analyse.get_sig_h('a') * c << endl;
    out << "Sigma (electron): " << analyse.get_sig_e('a') * c << endl;
    out << "Covariance (rh,re): " << analyse.get_cov('a') * c * c << endl;
    out << "Correlation factor (rh,re): " << analyse.get_corr('a') << endl;

    if (!aeqb) {

        out << "\n";
        out << "For the beta part:" << endl;
        out << "<rhx>,<rhy>,<rhz>: " << analyse.get_rh('x', 'b') * c << ", "
                << analyse.get_rh('y', 'b') * c << ", "
                << analyse.get_rh('z', 'b') * c << endl;
        out << "<rex>,<rey>,<rez>: " << analyse.get_re('x', 'b') * c << ", "
                << analyse.get_re('y', 'b') * c << ", "
                << analyse.get_re('z', 'b') * c << endl;
        out << "Electron hole separation |<rh-re>|: "
                << analyse.get_sep('b') * c << endl;
        out << "Distance for each coord: x= " << analyse.get_dex_c('x', 'b') * c
                << " y= " << analyse.get_dex_c('y', 'b') * c << " z="
                << analyse.get_dex_c('z', 'b') * c << endl;
        out << "Averaged distance over all coordinates: "
                << analyse.get_dex_tot('b') * c << endl;
        out << "Sigma (hole): " << analyse.get_sig_h('b') * c << endl;
        out << "Sigma (electron): " << analyse.get_sig_e('b') * c << endl;
        out << "Covariance (rh,re): " << analyse.get_cov('b') * c * c << endl;
        out << "Correlation factor (rh,re): " << analyse.get_corr('b') << endl;

    } else {

        out << "Alpha is equal to beta." << endl;

    }
    out << "\n";
    out <<"--------------------------------------------------------------------"
            << endl;
    out << endl;

}//end fct

}//end namespace libwfa
