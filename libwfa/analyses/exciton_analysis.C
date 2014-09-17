#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_analysis.h"

namespace libwfa {

using namespace arma;


exciton_analysis::exciton_analysis(const mom_builder_i &bld,
    const ab_matrix &tdm, size_t maxmm) :
    exciton_analysis_base(tdm.is_alpha_eq_beta(),
            std::min(bld.max_moment(), maxmm)) {

    calculate(bld, tdm.alpha(), moment(false));
    if (!tdm.is_alpha_eq_beta()) {
        calculate(bld, tdm.beta(), moment(true));
    }
}


void exciton_analysis::print_header(std::ostream &out, size_t off) const {

    std::string os(off, ' ');
    out << os << "Spatial multipole analysis of the transition density matrix"
            << std::endl;
    out << os << "in terms of the hole (r_h) and electron (r_e) coordinates."
            << std::endl;
}


void exciton_analysis::calculate(const mom_builder_i &bld,
    const arma::mat &tdm, exciton_moments &mom) {

    double n = bld.perform(tdm, 0, 0, 0, 0);

    for (size_t i = 1; i <= mom.n_max(); i++)
    for (size_t j = 0; j <= i; j++) {
        vec mj(3, fill::zeros);
        for (size_t k = 0; k < 3; k++)
            mj(k) = bld.perform(tdm, k, j, k, i - j) / n;
        mom.set(j, i - j, mj);
    }
}


void exciton_analysis::analysis(std::ostream &out,
        const exciton_moments &mom, size_t off) const {

    std::string os(off, ' ');
    out << std::setprecision(6) << std::fixed;
    { // Scope of linear quantities
        vec rh = mom.get(0, 1) * constants::au2ang;
        vec re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);
        out << os << "<r_h> [Ang]:" << std::string(24, ' ');
        print(out, rh);
        out << std::endl;
        out << os << "<r_e> [Ang]:" << std::string(24, ' ');
        print(out, re);
        out << std::endl;
        out << os << "|<r_e - r_h>| [Ang]:" << std::string(16, ' ')
                << std::setw(10) << tot << std::endl << std::endl;
    } // End of scope of rh, re, and tot

    { // Scope of quadratic quantities
        vec sh2 = mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1);
        vec se2 = mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0);
        vec d2  = mom.get(2, 0) + mom.get(0, 2) - mom.get(1, 1) * 2.;
        sh2 *= constants::au2ang * constants::au2ang;
        se2 *= constants::au2ang * constants::au2ang;
        d2  *= constants::au2ang * constants::au2ang;
        double sh = sqrt(accu(sh2)), se = sqrt(accu(se2));
        double cov = accu(mom.get(1, 1) - mom.get(1, 0) % mom.get(0, 1));
        cov *= constants::au2ang * constants::au2ang;

        //out << "<r_h,i^2> - <r_h,i>^2 [Ang^2]: ";
        //print(out, sh2);
        out << os << "Hole size [Ang]:" << std::string(20, ' ')
                << std::setw(10) << sh << std::endl;
        out << os << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(sh2));
        out << std::endl;
        //out << "<r_e,i^2> - <r_e,i>^2 [Ang^2]: ";
        //print(out, se2);
        out << os << "Electron size [Ang]:" << std::string(16, ' ')
                << std::setw(10) << se << std::endl;
        out << os << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(se2));
        out << std::endl;
        //out << "  <(r_e,i-r_h,i)^2> [Ang^2]:" << std::string(owidth-62, ' ');
        //print_vec(d2, out);
        out << os << "RMS electron-hole separation [Ang]: "
                << std::setw(10) << sqrt(accu(d2)) << std::endl;
        out << os << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(d2));
        out << std::endl;
        out << os << "Covariance(r_h, r_e) [Ang^2]:" << std::string(7, ' ')
                << std::setw(10) << cov << std::endl;
        out << os << "Correlation coefficient:" << std::string(12, ' ')
                << std::setw(10) << cov / (se * sh) << std::endl;
    }
}


}// end namespace libwfa

