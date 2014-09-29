#include <iomanip>
#include "ndo_analysis.h"

namespace libwfa {

using namespace arma;


ndo_analysis::ndo_analysis(const mat &s,
    const ab_matrix &c, const ab_matrix &ddm) {

    if (ddm.is_alpha_eq_beta()) {
        m_ndo[0] = new orbital_data(s, c.alpha(), ddm.alpha());
        m_ndo[1] = 0;
    }
    else {
        m_ndo[0] = new orbital_data(s, c.alpha(), ddm.alpha());
        m_ndo[1] = new orbital_data(s, c.beta(), ddm.beta());
    }
}


ndo_analysis::~ndo_analysis() {

    delete m_ndo[0]; m_ndo[0] = 0;
    if (m_ndo[1]) { delete m_ndo[1]; m_ndo[1] = 0; }
}


void ndo_analysis::form_ad(ab_matrix &at, ab_matrix &de) const {

    if (m_ndo[1]) {
        at.set_alpha_neq_beta();
        de.set_alpha_neq_beta();
        form_ad(m_ndo[0]->get_occ(), m_ndo[0]->get_coeff(),
                at.alpha(), de.alpha());
        form_ad(m_ndo[1]->get_occ(), m_ndo[1]->get_coeff(),
                at.beta(), de.beta());
    }
    else {
        at.set_alpha_eq_beta();
        de.set_alpha_eq_beta();
        form_ad(m_ndo[0]->get_occ(), m_ndo[0]->get_coeff(),
                at.alpha(), de.alpha());
    }
}


void ndo_analysis::analyse(std::ostream &out, size_t nndo) const {

    if (m_ndo[1]) {
        out << "NDOs (alpha):" << std::endl;
        analysis(out, m_ndo[0]->get_occ(), nndo);
        out << "NDOs (beta):" << std::endl;
        analysis(out, m_ndo[1]->get_occ(), nndo);
    }
    else {
        out << "NDOs:" << std::endl;
        analysis(out, m_ndo[0]->get_occ() * 2., nndo);
    }
}


void ndo_analysis::export_orbitals(orbital_printer_i &pr, double thresh) const {

    if (m_ndo[1]) {

        orbital_selector s_a, s_b;
        bld_selector(m_ndo[0]->get_occ(), thresh, s_a);
        bld_selector(m_ndo[1]->get_occ(), thresh, s_b);

        pr.perform(orbital_type::ndo, *m_ndo[0], s_a, *m_ndo[1], s_b);
    }
    else {

        orbital_selector s;
        bld_selector(m_ndo[0]->get_occ(), thresh, s);

        pr.perform(orbital_type::ndo, *m_ndo[0], s);
    }
}


void ndo_analysis::form_ad(const vec &e, const mat &c, mat &at, mat &de) {

    Col<uword> ix = find(e > 0.0, 1);
    if (ix.n_rows != 0) {
        if (ix(0) != 0) {
            mat ux = c.cols(ix(0), e.n_rows - 1);
            at = ux * diagmat(e.rows(ix(0), e.n_rows - 1)) * ux.t();
            ux = c.cols(0, ix(0) - 1);
            de = ux * diagmat(e.rows(0, ix(0) - 1)) * ux.t();
        }
        else {
            at = c * diagmat(e) * c.t();
            de = mat(c.n_cols, c.n_cols, fill::zeros);
        }
    }
    else {
        at = mat(c.n_cols, c.n_cols, fill::zeros);
        de = c * diagmat(e) * c.t();
    }
}


void ndo_analysis::analysis(std::ostream &out, const vec &ev, size_t nndo) {

    // Compute # attached and detached electrons first
    double na = 0.0, na2 = 0.0, nd = 0.0, nd2 = 0.0;

    size_t i = 0, nndo0 = 0;
    for (; i < ev.n_elem && ev(i) < 0; i++, nndo0++) {
        double cur = ev(i);
        nd += cur; nd2 += cur * cur;
    }
    nd *= -1;
    for (; i < ev.n_elem; i++) {
        double cur = ev(i);
        na += cur; na2 += cur * cur;
    }
    nndo0 = std::min(nndo0, ev.n_elem - nndo0);
    nndo = std::min(nndo, nndo0);

    std::string offset(2, ' ');
    out << "  Leading detachment eigenvalues: ";
    out << std::setprecision(4) << std::fixed;
    for (i = 0; i < nndo; i++) out << std::setw(9) << ev(i);
    out << std::endl;

    out << "  Leading attachment eigenvalues: ";
    out << std::setprecision(4) << std::fixed;
    size_t j = ev.n_elem - 1;
    for (i = 0; i < nndo; i++, j--) out << std::setw(9) << ev(j);
    out << std::endl;

    out << "  Number of detached / attached electrons: p_D = ";
    out << std::setw(7) << nd;
    out << ", p_A = " << std::setw(7) << na << std::endl;
    out << std::setprecision(6);
    out << "  Number of involved orbitals: PR_D = ";
    out << std::setw(9) << (nd * nd) / nd2;
    out << ", PR_A = ";
    out << std::setw(9) << (na * na) / na2 << std::endl;
}


void ndo_analysis::bld_selector(const arma::vec &e, double thresh,
    orbital_selector &sel) {

    size_t ntot = e.size();
    uvec p0 = find(e > -1. * thresh, 1), p1 = find(e > thresh, 1);

    size_t n0 = (size_t) (p0.size() == 1 ? p0(0) : 0);
    size_t n1 = (size_t) (p1.size() == 1 ? p1(0) : ntot);

    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    if (n0 > 0)    sel.select(true,   0,   n0, 1, true);
    if (n1 < ntot) sel.select(false, n1, ntot, 1, true);

}

} // namespace libwfa
