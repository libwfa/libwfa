#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_analysis_ad.h"

namespace libwfa {

using namespace arma;


exciton_analysis_ad::exciton_analysis_ad(const mom_builder_i &bld,
    const ab_matrix &adm, const ab_matrix &ddm, size_t maxmm) :
    exciton_analysis_base(adm.is_alpha_eq_beta(),
            std::min(bld.max_moment(), maxmm)) {

    calculate(bld, adm.alpha(), ddm.alpha(), moment(false));
    if (! adm.is_alpha_eq_beta()) {
        calculate(bld, adm.beta(), ddm.beta(), moment(true));
    }
}


void exciton_analysis_ad::print_header(std::ostream &out, size_t off) const {

    out << std::string(off, ' ');
    out << "Exciton analysis of the difference density matrix" << std::endl;
}


void exciton_analysis_ad::analysis(std::ostream &out,
    const exciton_moments &mom, size_t off) const {

    std::string os(off, ' ');
    out << std::setprecision(6) << std::fixed;
    { // Scope of rh, re, and tot
        vec rh = mom.get(0, 1) * constants::au2ang;
        vec re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);

        out << os << "<r_h> [Ang]:" << std::string(18, ' ');
        print(out, rh);
        out << std::endl;
        out << os << "<r_e> [Ang]:" << std::string(18, ' ');
        print(out, re);
        out << std::endl;
        out << os << "|<r_e - r_h>| [Ang]:" << std::string(10, ' ')
                << std::setw(10) << tot << std::endl;
    } // End of scope of rh, re, and tot

    { // Scope of quadratic quantities
        vec sh2 = mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1);
        vec se2 = mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0);
        sh2 *= constants::au2ang * constants::au2ang;
        se2 *= constants::au2ang * constants::au2ang;
        double sh = sqrt(accu(sh2)), se = sqrt(accu(se2));

        //out << "  <r_h,i^2> - <r_h,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(sh2, out);
        out << os << "Hole size [Ang]:" << std::string(14, ' ')
            << std::setw(10) << sh << std::endl;
        out << os << "  Cartesian components [Ang]: ";
        print(out, sqrt(sh2));
        out << std::endl;
        //out << "  <r_e,i^2> - <r_e,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(se2, out);
        out << os << "Electron size [Ang]:" << std::string(10, ' ')
            << std::setw(10) << se << std::endl;
        out << os << "  Cartesian components [Ang]: ";
        print(out, sqrt(se2));
        out << std::endl;
    }
}


void exciton_analysis_ad::calculate(const mom_builder_i &bld,
    const arma::mat &adm, const arma::mat &ddm, exciton_moments &mom) {

    double n = bld.perform(adm, 'x', 0);
    vec m0(3, fill::zeros);
    m0(0) = n;
    mom.set(0, 0, m0);

    for (size_t i = 1; i <= mom.n_max(); i++) {
        vec mh(3, fill::zeros), me(3, fill::zeros);
        for (size_t k = 0; k < 3; k++) {
            me(k) = bld.perform(adm, k, i) / n;
            mh(k) = bld.perform(ddm * -1., k, i) / n;
        }
        mom.set(i, 0, me);
        mom.set(0, i, mh);
    }
}


}// end namespace libwfa

