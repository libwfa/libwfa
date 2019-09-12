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
    if (accu(abs(mom.get(0, 0)))<1.e-6) {
        out << os << "... vanishing." << std::endl;
        return;
    }
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

void exciton_analysis_ad::combine(const exciton_moments &a,
    const exciton_moments &b, exciton_moments &res) const {

    size_t nmax = std::min(a.n_max(), b.n_max());
    if (nmax != res.n_max()) {
        throw libwfa_exception("exciton_analysis_base", "combine",
                __FILE__, __LINE__, "nmax");
    }

    // Normalization factors: alpha/beta and attachment/detachment
    vec m0a = a.get(0,0);
    vec m0b = b.get(0,0);

    double paa = m0a(1), pad = m0a(2); //temp!
    double pba = m0b(1), pbd = m0b(2);

    if (paa + pba == 0.0) {
        throw libwfa_exception("exciton_analysis_base", "combine",
                __FILE__, __LINE__, "paa+pba==0, attachment");
    }

    if (pad + pbd == 0.0) {
        throw libwfa_exception("exciton_analysis_base", "combine",
                __FILE__, __LINE__, "pad+pbd==0, detachment");
    }

    for (size_t i = 0; i <= nmax; i++) {
        { // attachment
            vec maa = a.get(i, 0), mba = b.get(i, 0);
            vec mc = (paa * maa + pba * mba) / (paa + pba);
            res.set(i, 0, mc);
        }

        { // detachment
            vec mad = a.get(0, i), mbd = b.get(0, i);
            vec mc = (pad * mad + pbd * mbd) / (pad + pbd);
            res.set(0, i, mc);
        }
    }
}

void exciton_analysis_ad::calculate(const mom_builder_i &bld,
    const arma::mat &adm, const arma::mat &ddm, exciton_moments &mom) {

    double pa = bld.perform(adm, 'x', 0);
    double pd = bld.perform(ddm, 'x', 0);
    vec m0(3, fill::zeros);
    m0(1) = pa;
    m0(2) = pd;
    mom.set(0, 0, m0);

    for (size_t i = 1; i <= mom.n_max(); i++) {
        vec md(3, fill::zeros), ma(3, fill::zeros);
        for (size_t k = 0; k < 3; k++) {
            ma(k) = bld.perform(adm, k, i) / pa;
            md(k) = bld.perform(ddm, k, i) / pd;
        }
        mom.set(i, 0, ma);
        mom.set(0, i, md);
    }
}


}// end namespace libwfa

