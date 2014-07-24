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

    const size_t owidth = 80;
    
    out << " +" << std::string(owidth, '-') << "+" << std::endl;
    out << " |" << std::string(10, ' ') 
        << "Spatial multipole analysis of the difference density matrix"
        << std::string(11, ' ') << "|" << std::endl
        << " |" << std::string(10, ' ')
        << " in terms of the hole (r_h) and electron (r_e) coordinates."
        << std::string(11, ' ') << "|" << std::endl;
    out << " +" << std::string(owidth, '-') << "+" << std::endl;    
    
    out << std::setprecision(6) << std::fixed;
    { // Scope of rh, re, and tot
        Col<double> rh = mom.get(0, 1) * constants::au2ang;
        Col<double> re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);
                
        out << "  <r_h> [Ang]:" << std::string(owidth-48, ' ');
        print_vec(rh, out);

        out << "  <r_e> [Ang]:" << std::string(owidth-48, ' ');
        print_vec(re, out);
        
        out << "  |<r_e - r_h>| [Ang]:" << std::string(owidth-30, ' ')
            << std::setw(10) << tot << std::endl << std::endl;        
    } // End of scope of rh, re, and tot

    { // Scope of quadratic quantities
        Col<double> sh2 = mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1);
        sh2 *= constants::au2ang * constants::au2ang;
        double sh = sqrt(accu(sh2));
        
        out << "  <r_h,i^2> - <r_h,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        print_vec(sh2, out);
        out << "  Hole size [Ang]:"  << std::string(owidth-26, ' ')
            << std::setw(10) << sh << std::endl;
            
        Col<double> se2 = mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0);
        se2 *= constants::au2ang * constants::au2ang;
        double se = sqrt(accu(se2));
            
        out << "  <r_e,i^2> - <r_e,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        print_vec(se2, out);
        out << "  Electron size [Ang]:"  << std::string(owidth-30, ' ')
            << std::setw(10) << se << std::endl;        
    }    
    out << " +" << std::string(owidth, '-') << "+" << std::endl;    
}


}//end namespace
