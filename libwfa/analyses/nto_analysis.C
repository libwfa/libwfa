#include <iomanip>
#include "nto_analysis.h"

namespace libwfa {

using namespace arma;


nto_analysis::nto_analysis(const mat &s, const ab_matrix &c,
    const ab_matrix &tdm) {

    ab_matrix edm, hdm;
    form_eh(s, tdm, edm, hdm);

    initialize(s, c, edm, hdm);
}


nto_analysis::nto_analysis(const mat &s, const ab_matrix &c,
    const ab_matrix &edm, const ab_matrix &hdm) {

    initialize(s, c, edm, hdm);
}


void nto_analysis::analyse(std::ostream &out, size_t nnto) const {

    if (nnto == 0)  return;
    if (m_nto[2]) {
        out << "NTOs (alpha-electron)" << std::endl;
        analysis(out, m_nto[0]->get_occ(), nnto);
        out << "NTOs (alpha-hole)" << std::endl;
        analysis(out, m_nto[1]->get_occ(), nnto);
        out << "NTOs (beta-electron)" << std::endl;
        analysis(out, m_nto[0]->get_occ(), nnto);
        out << "NTOs (beta-hole)" << std::endl;
        analysis(out, m_nto[1]->get_occ(), nnto);
    }
    else {
        out << "NTOs (electron)" << std::endl;
        analysis(out, m_nto[0]->get_occ() * 2.0, nnto);
        out << "NTOs (hole)" << std::endl;
        analysis(out, m_nto[1]->get_occ() * 2.0, nnto);
    }
}


void nto_analysis::export_orbitals(orbital_printer_i &pr, double thresh) const {

    if (m_nto[2]) {

        const vec &ee_a = m_nto[0]->get_occ(),   &eh_a = m_nto[1]->get_occ();
        const mat &ce_a = m_nto[0]->get_coeff(), &ch_a = m_nto[1]->get_coeff();
        const vec &ee_b = m_nto[2]->get_occ(),   &eh_b = m_nto[2]->get_occ();
        const mat &ce_b = m_nto[3]->get_coeff(), &ch_b = m_nto[3]->get_coeff();
        orbital_data nto_a(join_cols(eh_a * -1., flipud(ee_a)),
                join_rows(ch_a, fliplr(ce_a)));
        orbital_data nto_b(join_cols(eh_b * -1., flipud(ee_b)),
                join_rows(ch_b, fliplr(ce_b)));
        orbital_selector s_a, s_b;
        build_selector(ee_a, eh_a, thresh, s_a);
        build_selector(ee_b, eh_b, thresh, s_b);

        pr.perform(orbital_type::nto, nto_a, s_a, nto_b, s_b);
    }
    else {

        const vec &ee = m_nto[0]->get_occ(), &eh = m_nto[1]->get_occ();
        const mat &ce = m_nto[0]->get_coeff(), &ch = m_nto[1]->get_coeff();
        orbital_data nto(join_cols(eh * -1., flipud(ee)),
                join_rows(ch, fliplr(ce)));
        orbital_selector s;
        build_selector(ee, eh, thresh, s);

        pr.perform(orbital_type::nto, nto, s);
    }
}


void nto_analysis::form_eh(const mat &s, const ab_matrix &tdm,
    ab_matrix &edm, ab_matrix &hdm) {

    if (tdm.is_alpha_eq_beta()) {
        edm.set_alpha_eq_beta();
        hdm.set_alpha_eq_beta();

        form_eh(s, tdm.alpha(), edm.alpha(), hdm.alpha());
    }
    else {
        edm.set_alpha_neq_beta();
        hdm.set_alpha_neq_beta();

        form_eh(s, tdm.alpha(), edm.alpha(), hdm.alpha());
        form_eh(s, tdm.beta(), edm.beta(), hdm.beta());
    }
}


void nto_analysis::initialize(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm) {

    if (edm.is_alpha_eq_beta()) {
        m_nto[0] = new orbital_data(s, c.alpha(), edm.alpha());
        m_nto[1] = new orbital_data(s, c.alpha(), hdm.alpha());
        m_nto[2] = m_nto[3] = 0;
    }
    else {
        m_nto[0] = new orbital_data(s, c.alpha(), edm.alpha());
        m_nto[1] = new orbital_data(s, c.alpha(), hdm.alpha());
        m_nto[2] = new orbital_data(s, c.beta(), edm.beta());
        m_nto[3] = new orbital_data(s, c.beta(), hdm.beta());
    }
}


void nto_analysis::analysis(std::ostream &out,
    const arma::vec &e, size_t nnto) {

    out << "  Leading SVs:" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "  ";
    for (size_t i = 0, j = e.n_rows - 1; i < nnto; i++, j--)
        out << std::setw(9) << e(j);
    out << std::endl;

    double total = accu(e);
    out << std::setprecision(6) << std::fixed;
    out << "  Sum of SVs:  " << std::setw(11) << total << std::endl;
    out << "  Participation ratio (PR_NTO):  "
        << std::setw(11) << total * total / dot(e, e);
    out << std::endl;
}


void nto_analysis::build_selector(const arma::vec &e, const arma::vec &h,
    double thresh, orbital_selector &sel) {

    size_t ntot = h.size() + e.size();
    uvec ph = find(h > thresh, 1), pe = find(e > thresh, 1);

    size_t nh = (ph.size() == 1 ? ph(0) : 0);
    size_t ne = (pe.size() == 1 ? pe(0) : 0);
    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    if (nh > 0) sel.select(true, nh, h.size(), 1);
    if (ne > 0) sel.select(false, h.size(), h.size() + (e.size() - ne), 1);
}


} // namespace libwfa


