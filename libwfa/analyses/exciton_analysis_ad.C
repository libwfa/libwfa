#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_analysis_ad.h"

namespace libwfa {

using namespace arma;


exciton_analysis_ad::exciton_analysis_ad(const mom_builder_i &bld,
    const ab_matrix &adm, const ab_matrix &ddm, size_t maxmm) {

    maxmm = std::min(bld.max_moment(), maxmm);

    m_mom[0] = new exciton_moments(maxmm);
    calculate(bld, adm.alpha(), ddm.alpha(), *m_mom[0]);

    if (adm.is_alpha_eq_beta()) { m_mom[1] = 0; }
    else {
        m_mom[1] = new exciton_moments(maxmm);
        calculate(bld, adm.beta(), ddm.beta(), *m_mom[1]);
    }
}


exciton_analysis_ad::~exciton_analysis_ad() {

    delete m_mom[0];
    if (m_mom[1]) delete m_mom[1];
}


void exciton_analysis_ad::analyse(std::ostream &out) const {

    if (m_mom[1]) {
        out << "alpha spin:" << std::endl;
        analysis(out, *m_mom[0]);
        out << "beta spin:" << std::endl;
        analysis(out, *m_mom[1]);
    }
    else {
        analysis(out, *m_mom[0]);
    }
}


void exciton_analysis_ad::calculate(const mom_builder_i &bld,
    const arma::mat &adm, const arma::mat &ddm, exciton_moments &mom) {

    double n = bld.perform(adm, 'x', 0);

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


void exciton_analysis_ad::analysis(std::ostream &out,
    const exciton_moments &mom) {

    out << std::setprecision(6) << std::fixed;
    { // Scope of rh, re, and tot
        vec rh = mom.get(0, 1) * constants::au2ang;
        vec re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);

        out << "<r_h> [Ang]:" << std::string(18, ' ');
        print(out, rh);
        out << std::endl;
        out << "<r_e> [Ang]:" << std::string(18, ' ');
        print(out, re);
        out << std::endl;
        out << "|<r_e - r_h>| [Ang]:" << std::string(10, ' ')
                << std::setw(10) << tot << std::endl << std::endl;
    } // End of scope of rh, re, and tot

    { // Scope of quadratic quantities
        vec sh2 = mom.get(0, 2) - mom.get(0, 1) % mom.get(0, 1);
        vec se2 = mom.get(2, 0) - mom.get(1, 0) % mom.get(1, 0);
        sh2 *= constants::au2ang * constants::au2ang;
        se2 *= constants::au2ang * constants::au2ang;
        double sh = sqrt(accu(sh2)), se = sqrt(accu(se2));

        //out << "  <r_h,i^2> - <r_h,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(sh2, out);
        out << "Hole size [Ang]:" << std::string(10, ' ')
            << std::setw(10) << sh << std::endl;
        out << "  Cartesian components [Ang]: ";
        print(out, sqrt(sh2));
        out << std::endl;
        //out << "  <r_e,i^2> - <r_e,i>^2 [Ang^2]:" << std::string(owidth-66, ' ');
        //print_vec(se2, out);
        out << "Electron size [Ang]:" << std::string(10, '-')
            << std::setw(10) << se << std::endl;
        out << "  Cartesian components [Ang]: ";
        print(out, sqrt(se2));
        out << std::endl;
    }
}


void exciton_analysis_ad::print(std::ostream &out,
    const vec &vec, size_t width) {

    out << "[";
    if (vec.n_rows > 0) out << std::setw(width) << vec(0);
    for (size_t i = 1; i < vec.n_rows; i++) {
        out << ", " << std::setw(width) << vec(i);
    }
    out << "]";
}


}// end namespace libwfa





