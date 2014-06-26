#include <iomanip>
#include <fstream>
#include "ex_ana_printer_ad.h"

namespace libwfa{

using namespace arma;

const double ex_ana_printer_ad::k_au2ang = 0.52917724900;

void ex_ana_printer_ad::perform(ex_analyse_ad &analyse_ad, std::ostream &out){

    out << std::setw(10) << std::setprecision(6) << std::fixed;
    out << std::endl << std::endl;
    out << std::string(70, '-') << std::endl;
    out << "Performing exciton analysis (Att//Det)... (in Ang) " << std::endl;
    out << std::endl;
    out << "For the alpha part:" << std::endl;
    do_print('a', analyse_ad, out);
    if (! analyse_ad.aeqb()) {
        out << std::endl;
        out << "For the beta part:" << std::endl;
        do_print('b', analyse_ad, out);
    } else {
        out << "Alpha is equal to beta." << std::endl;
    }
    out << std::endl;
    out << std::string(70, '-') << std::endl;
    out << std::endl;

}//end fct

void ex_ana_printer_ad::do_print(char spin, ex_analyse_ad &analyse_ad, std::ostream &out){

    out << "<rhx>,<rhy>,<rhz>: "
            << analyse_ad.get_rh('x', spin) * k_au2ang << ", "
            << analyse_ad.get_rh('y', spin) * k_au2ang << ", "
            << analyse_ad.get_rh('z', spin) * k_au2ang
            << std::endl;
    out << "<rex>,<rey>,<rez>: "
            << analyse_ad.get_re('x', spin) * k_au2ang << ", "
            << analyse_ad.get_re('y', spin) * k_au2ang << ", "
            << analyse_ad.get_re('z', spin) * k_au2ang
            << std::endl;
    out << "Electron hole separation |<rh-re>|: "
            << analyse_ad.get_sep(spin) * k_au2ang
            << std::endl;
    out << "Sigma (hole): "
            << analyse_ad.get_sig_h(spin) * k_au2ang << std::endl;
    out << "Sigma (electron): "
            << analyse_ad.get_sig_e(spin) * k_au2ang << std::endl;


}//end fct

}//end namespace
