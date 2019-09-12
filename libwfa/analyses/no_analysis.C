#include <iomanip>
#include "no_analysis.h"

namespace libwfa {

using namespace arma;


no_analysis::no_analysis(const mat &s, const ab_matrix &c,
    const ab_matrix &sdm) {

    m_no[0] = new orbital_data(s, c.alpha(), sdm.alpha());
    if (sdm.is_alpha_eq_beta()) {
        m_no[1] = m_no[2] = 0;
    }
    else {
        m_no[1] = new orbital_data(s, c.alpha(), sdm.beta());
        m_no[2] = new orbital_data(s, c.alpha(), sdm.alpha() + sdm.beta());
    }
}


no_analysis::~no_analysis() {

    delete m_no[0]; m_no[0] = 0;
    if (m_no[1]) { delete m_no[1]; m_no[1] = 0; }
    if (m_no[2]) { delete m_no[2]; m_no[2] = 0; }
}


void no_analysis::analyse(std::ostream &out, size_t nno) const {

    if (nno == 0)  return;
    if (m_no[1]) {
        out << "NOs (alpha)" << std::endl;
        analysis_p1(out, m_no[0]->get_occ(), nno);
        out << "NOs (beta)" << std::endl;
        analysis_p1(out, m_no[1]->get_occ(), nno);

        out << "NOs (spin-traced)" << std::endl;
        analysis_p1(out, m_no[2]->get_occ(), nno);
        analysis_p2(out, m_no[2]->get_occ());
    }
    else {
        vec ev(m_no[0]->get_occ());
        ev *= 2.0;
        out << "NOs" << std::endl;
        analysis_p1(out, ev, nno);
        analysis_p2(out, ev);
    }
}


void no_analysis::export_orbitals(orbital_printer_i &pr, double thresh) const {

    if (m_no[1]) {
        orbital_selector s_a, s_b;
        build_selector(m_no[0]->get_occ(), thresh, s_a);
        build_selector(m_no[1]->get_occ(), thresh, s_b);
        pr.perform(orbital_type::no, *m_no[0], s_a, *m_no[1], s_b);
    }
    else {
        orbital_selector s;
        build_selector(m_no[0]->get_occ(), thresh, s);
        pr.perform(orbital_type::no, *m_no[0], s);
    }
}


void no_analysis::analysis_p1(std::ostream &out, const vec &ev, size_t nno) {

    double nelec = accu(ev);
    size_t ihomo = ev.n_elem - (size_t)(nelec / ev(ev.n_elem - 1) + 0.5);

    size_t min = (ihomo > nno ? ihomo - nno : 0);
    size_t max = (ihomo + nno < ev.n_elem ? ihomo + nno : ev.n_elem);

    out << "  Occupation of frontier NOs:" << std::endl;
    out << std::fixed << std::setprecision(4);
    out << "  ";
    for (size_t i = min; i < max; i++) out << std::setw(9) << ev(i);
    out << std::endl;
    out << "  Number of electrons: ";
    out << std::setprecision(6) << std::setw(9) << nelec << std::endl;
}


void no_analysis::analysis_p2(std::ostream &out, const vec &ev) {

    double nu = 0.0, nu2 = 0.0, nunl = 0.0;
    for (size_t i = 0; i < ev.n_elem; i++) {

        double n = ev(i), nn = 2. - n;
        nunl += n * n * nn * nn;

        n = std::min(n, nn);
        nu += n;
        nu2 += n * n;
    }

    out << std::setprecision(5);
    out << "  Number of unpaired electrons: n_u = ";
    out << std::setw(8) << nu;
    out << ", n_u,nl = ";
    out << std::setw(8) << nunl << std::endl;

    if (nu2 > 1.e-6) {
        out << "  NO participation ratio (PR_NO): ";
        out << std::setw(9) << std::setprecision(6) << (nu * nu) / nu2 << std::endl;
    }
}


void no_analysis::build_selector(const vec &e, double thresh,
    orbital_selector &sel) {

    size_t ntot = e.size(), nelec = accu(e), pos = 0;
    uvec v = find(e > thresh, 1);
    if (v.size() != 0) pos = v(0);

    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    size_t nv = std::max(pos, ntot - 2 * nelec);
    size_t no = std::max(pos, ntot - nelec);
    sel.select(true, no, ntot, 1, true);
    if (nv < no) sel.select(false, nv, no, 1, true);
}

} // namespace libwfa


