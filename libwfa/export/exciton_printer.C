#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_printer.h"

namespace libwfa{

using namespace arma;

void exciton_printer::perform(ab_exciton_moments &mom,
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

void exciton_printer::print(exciton_moments &mom, std::ostream &out) const {
    
    const size_t owidth = 80;

    out << " +" << std::string(owidth, '-') << "+" << std::endl;
    out << " |" << std::string(10, ' ') 
        << "Spatial multipole analysis of the transition density matrix"
        << std::string(11, ' ') << "|" << std::endl
        << " |" << std::string(10, ' ')
        << " in terms of the hole (r_h) and electron (r_e) coordinates."
        << std::string(11, ' ') << "|" << std::endl;
    out << " +" << std::string(owidth, '-') << "+" << std::endl;
    
    out << std::setprecision(6) << std::fixed;
    { // Scope of linear quantities
        vec rh = mom.get(0, 1) * constants::au2ang;
        vec re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);
        out << "  <r_h> [Ang]:" << std::string(owidth-48, ' ');
        print_vec(rh, out);

        out << "  <r_e> [Ang]:" << std::string(owidth-48, ' ');
        print_vec(re, out);
        
        out << "  |<r_e - r_h>| [Ang]:" << std::string(owidth-30, ' ')
            << std::setw(10) << tot << std::endl << std::endl;
    } // End of scope of rh, re, and tot

    { // Scope of quadratic quantities
        vec sh2 = mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1);
        sh2 *= constants::au2ang * constants::au2ang;
        double sh = sqrt(accu(sh2));
        vec ssh2 = sqrt(sh2);
        
        //out << "  <r_h,i^2> - <r_h,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(sh2, out);
        out << "  Hole size [Ang]:"  << std::string(owidth-26, ' ')
            << std::setw(10) << sh << std::endl;
        out << "    Cartesian components [Ang]:" << std::string(owidth-65, ' ');
        print_vec (ssh2, out);
            
        vec se2 = mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0);
        se2 *= constants::au2ang * constants::au2ang;
        double se = sqrt(accu(se2));
        vec sse2 = sqrt(se2);
            
        //out << "  <r_e,i^2> - <r_e,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(se2, out);
        out << "  Electron size [Ang]:"  << std::string(owidth-30, ' ')
            << std::setw(10) << se << std::endl;
        out << "    Cartesian components [Ang]:" << std::string(owidth-65, ' ');
        print_vec (sse2, out);
        out << std::endl;
        
        vec d2 = mom.get(2,0) + mom.get(0,2) - mom.get(1,1) * 2.;
        d2 *= constants::au2ang * constants::au2ang;
        vec sd2 = sqrt(d2);        
        
        //out << "  <(r_e,i-r_h,i)^2> [Ang^2]:" << std::string(owidth-62, ' ');
        //print_vec(d2, out);
        out << "  RMS electron-hole separation [Ang]:"  << std::string(owidth-45, ' ')
            << std::setw(10) << sqrt(accu(d2)) << std::endl;
        out << "    Cartesian components [Ang]:" << std::string(owidth-65, ' ');
        print_vec (sd2, out);
                    
        double cov, cf;
        cov = accu(mom.get(1, 1) - mom.get(1, 0) % mom.get(0, 1));
        cov *= constants::au2ang * constants::au2ang;
        cf = cov / (se * sh);
        out << "  Covariance(r_h, r_e) [Ang^2]:" << std::string(owidth-39, ' ')
            << std::setw(10) << cov << std::endl;
        out << "  Correlation coefficient:" << std::string(owidth-34, ' ')
            << std::setw(10) << cf << std::endl;
    }
    out << " +" << std::string(owidth, '-') << "+" << std::endl;
}


}//end namespace libwfa
