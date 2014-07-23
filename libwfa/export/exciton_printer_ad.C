#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_printer_ad.h"

namespace libwfa{

using namespace arma;

void exciton_printer_ad::perform(ab_exciton_moments &mom,
    std::ostream &out) const {

    if (mom.is_alpha_eq_beta()) {
        print(mom.alpha(), out);
    }
    else {
        out << "alpha spin:" << std::endl;
        print(mom.alpha(), out);
        out << "beta spin:" << std::endl;
        print(mom.beta(), out);
    }
}


void exciton_printer_ad::print(exciton_moments &mom, std::ostream &out) const {

    out << std::setprecision(6) << std::fixed;
    { // Scope of rh, re, and tot
        Col<double> rh = mom.get(0, 1);
        Col<double> re = mom.get(1, 0);
        double tot = norm(rh - re);
        out << "Hole center (in Ang):     [";
        for (size_t i = 0; i < 3; i++) {
            out << " " << std::setw(10) << rh(i) * constants::au2ang;
            if (i != 2) out << ",";
        }
        out << "]" << std::endl;

        out << "Electron center (in Ang): [";
        for (size_t i = 0; i < 3; i++) {
            out << " " << std::setw(10) << re(i) * constants::au2ang;
            if (i != 2) out << ",";
        }
        out << "]" << std::endl;
        out << "Electron-hole distance:   ";
        out << std::setw(10) << tot * constants::au2ang << " Ang" << std::endl;
    } // End of scope of rh, re, and tot

    { // Scope of sh, se, cov, cf
        double sh, se;
        sh = sqrt(fabs(accu(mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1))));
        se = sqrt(fabs(accu(mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0))));
        out << "Variance of the hole:     ";
        out << std::setw(10) << sh * constants::au2ang << " Ang" << std::endl;
        out << "Variance of the electron: ";
        out << std::setw(10) << se * constants::au2ang << " Ang" << std::endl;
    } // End scope of sh, se
}


}//end namespace
