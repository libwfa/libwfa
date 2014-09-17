#include <iomanip>
#include <libwfa/core/constants.h>
#include "exciton_analysis.h"

namespace libwfa {

using namespace arma;


exciton_analysis::exciton_analysis(const mom_builder_i &bld,
    const ab_matrix &tdm, size_t maxmm) {

    maxmm = std::min(bld.max_moment(), maxmm);

    m_mom[0] = new exciton_moments(maxmm);
    calculate(bld, tdm.alpha(), *m_mom[0]);

    if (tdm.is_alpha_eq_beta()) { m_mom[1] = 0; }
    else {
        m_mom[1] = new exciton_moments(maxmm);
        calculate(bld, tdm.beta(),  *m_mom[1]);
    }
}


exciton_analysis::~exciton_analysis() {

    delete m_mom[0];
    if (m_mom[1]) delete m_mom[1];
}


void exciton_analysis::analyse(std::ostream &out) const {

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


void exciton_analysis::analysis(std::ostream &out, const exciton_moments &mom) {

    out << std::setprecision(6) << std::fixed;
    { // Scope of linear quantities
        vec rh = mom.get(0, 1) * constants::au2ang;
        vec re = mom.get(1, 0) * constants::au2ang;
        double tot = norm(rh - re);
        out << "<r_h> [Ang]:" << std::string(24, ' ');
        print(out, rh);
        out << std::endl;
        out << "<r_e> [Ang]:" << std::string(24, ' ');
        print(out, re);
        out << std::endl;
        out << "|<r_e - r_h>| [Ang]:" << std::string(16, ' ')
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
        out << "Hole size [Ang]:" << std::string(20, ' ')
                << std::setw(10) << sh << std::endl;
        out << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(sh2));
        out << std::endl;
        //out << "<r_e,i^2> - <r_e,i>^2 [Ang^2]: ";
        //print(out, se2);
        out << "Electron size [Ang]:" << std::string(16, ' ')
                << std::setw(10) << se << std::endl;
        out << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(se2));
        out << std::endl;
        //out << "  <(r_e,i-r_h,i)^2> [Ang^2]:" << std::string(owidth-62, ' ');
        //print_vec(d2, out);
        out << "RMS electron-hole separation [Ang]: "
                << std::setw(10) << sqrt(accu(d2)) << std::endl;
        out << "  Cartesian components [Ang]:" << std::string(7, ' ');
        print(out, sqrt(d2));
        out << std::endl;
        out << "Covariance(r_h, r_e) [Ang^2]:" << std::string(7, ' ')
                << std::setw(10) << cov << std::endl;
        out << "Correlation coefficient:" << std::string(12, ' ')
                << std::setw(10) << cov / (se * sh) << std::endl;
    }
}


void exciton_analysis::print(std::ostream &out, const vec &vec, size_t width) {

    out << "[";
    if (vec.n_rows > 0) out << std::setw(width) << vec(0);
    for (size_t i = 1; i < vec.n_rows; i++) {
        out << ", " << std::setw(width) << vec(i);
    }
    out << "]";
}


}// end namespace libwfa

