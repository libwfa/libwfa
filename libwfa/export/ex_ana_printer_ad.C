#include <iomanip>
#include <fstream>
#include "ex_ana_printer_ad.h"

namespace libwfa{

using namespace arma;

void ex_ana_printer_ad::perform(bool aeqb, ex_analyse_ad &analyse_ad, std::ostream &out){

    do_print (aeqb, analyse_ad, out);

}//end fct

void ex_ana_printer_ad::do_print(bool aeqb, ex_analyse_ad &analyse_ad, std::ostream &out){
    double c=0.52917724900;
    out<< std::setw(10) << std::setprecision(6) << std::fixed;
    out<<"\n \n";
    out<<"---------------------------------------------------------------------"
            << endl;
    out << "Performing exciton analysis (Attachement//Detachment)... "
            "(in \u212B) " << endl;

    out << "\n";
    out << "For the alpha part:" << endl;
    out << "<rhx>,<rhy>,<rhz>: " << analyse_ad.get_rh('x', 'a') * c << ", "
            << analyse_ad.get_rh('y', 'a') * c << ", "
            << analyse_ad.get_rh('z', 'a') * c << endl;
    out << "<rex>,<rey>,<rez>: " << analyse_ad.get_re('x', 'a') * c << ", "
            << analyse_ad.get_re('y', 'a') * c << ", "
            << analyse_ad.get_re('z', 'a') * c << endl;
    out << "Electron hole separation |<rh-re>|: " << analyse_ad.get_sep('a') * c
            << endl;
    out << "Sigma (hole): " << analyse_ad.get_sig_h('a') * c << endl;
    out << "Sigma (electron): " << analyse_ad.get_sig_e('a') * c << endl;


    if (!aeqb) {

        out << "\n";
        out << "For the beta part:" << endl;
        out << "<rhx>,<rhy>,<rhz>: " << analyse_ad.get_rh('x', 'b') * c << ", "
                << analyse_ad.get_rh('y', 'b') * c << ", "
                << analyse_ad.get_rh('z', 'b') * c << endl;
        out << "<rex>,<rey>,<rez>: " << analyse_ad.get_re('x', 'b') * c << ", "
                << analyse_ad.get_re('y', 'b') * c << ", "
                << analyse_ad.get_re('z', 'b') * c << endl;
        out << "Electron hole separation |<rh-re>|: "
                << analyse_ad.get_sep('b') * c << endl;
        out << "Sigma (hole): " << analyse_ad.get_sig_h('b') * c << endl;
        out << "Sigma (electron): " << analyse_ad.get_sig_e('b') * c << endl;


    } else {

        out << "Alpha is equal to beta." << endl;

    }
    out << "\n";
    out <<"--------------------------------------------------------------------"
            << endl;
    out << endl;

}//end fct

}//end namespace
