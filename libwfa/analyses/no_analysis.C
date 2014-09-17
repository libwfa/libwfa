#include <iomanip>
#include "no_analysis.h"

namespace libwfa {

using namespace arma;


/*
    std::ostream out;
    mat s;
    mat c;

    Steps to perform for the NO analysis given a density matrix as ab_matrix

    ab_matrix sdm;
    if (sdm.is_alpha_eq_beta()) {
        orbital_data no(s, c, sdm.alpha() * 2);
        analyse(no, out, 5);
        analyse_nu(no, out);
        export_orbitals(no);
    }
    else {
        orbital_data no_a(s, c, sdm.alpha());
        orbital_data no_b(s, c, sdm.beta());
        orbital_data no_x(s, c, sdm.alpha() + sdm.beta());
        out << "Alpha NOs" << std::endl;
        analyse(no_a, out, 5);
        out << "Beta NOs" << std::endl;
        analyse(no_b, out, 5);
        analyse_nu(no_x, out);
        export_orbitals(no_a, no_b);
    }
    export_density(sdm);
    pop_analysis(sdm);

    Steps to perform for the NDO analysis given a difference density matrix as ab_matrix

    ab_matrix ddm, atdm, dedm;
    if (ddm.is_alpha_eq_beta()) {
        orbital_data ndo(s, c, ddm.alpha() * 2);
        analyse(ndo, out, 5);
        export_orbitals(ndo);
        form_ad(ndo, atdm.alpha(), dedm.alpha());
    }
    else {
        orbital_data ndo_a(s, c, ddm.alpha());
        orbital_data ndo_b(s, c, ddm.beta());
        out << "Alpha NDOs" << std::endl;
        analyse(ndo_a, out, 5);
        out << "Beta NDOs" << std::endl;
        analyse(ndo_b, out, 5);
        export_orbitals(ndo_a, ndo_b);
        form_ad(ndo_a, atdm.alpha(), dedm.alpha());
        form_ad(ndo_b, atdm.beta(), dedm.beta());
    }
    export_density(ddm);
    export_density(atdm);
    export_density(dedm);
    pop_analysis(atdm, dedm);

    Steps to perform for the NTO analysis given a transition density matrix as ab_matrix

    ab_matrix tdm, edm, hdm;
    ctnumbers(s, tdm);

    form_eh(s, tdm, edm, hdm);
    if (tdm.is_alpha_eq_beta()) {
        orbital_data nto_e(s, c, edm.alpha() * 2);
        orbital_data nto_h(s, c, hdm.alpha() * 2);
        analyse(nto_e, out, 5);
        analyse(nto_h, out, 5);
        export_orbitals(nto_e + nto_h);
    }
    else {
        orbital_data nto_ea(s, c, edm.alpha() * 2);
        orbital_data nto_eb(s, c, edm.beta() * 2);
        orbital_data nto_ha(s, c, hdm.alpha() * 2);
        orbital_data nto_hb(s, c, hdm.beta() * 2);
        analyse(nto_ea, out, 5);
        analyse(nto_ha, out, 5);
        analyse(nto_eb, out, 5);
        analyse(nto_hb, out, 5);
        export_orbitals(nto_ea + nto_ha, nto_eb + nto_hb);
    }
    export_density(tdm);
    export_density(edm);
    export_density(hdm);

    Steps to perform for the SA-NTO analysis given a transition density matrix as ab_matrix

    ab_matrix tdm, edm_av, hdm_av;
    if (tdm.is_alpha_eq_beta()) {
        orbital_data nto_e(s, c, edm_av.alpha() * 2);
        orbital_data nto_h(s, c, hdm_av.alpha() * 2);
        analyse(nto_e, out, 5);
        analyse(nto_h, out, 5);
        export_orbitals(nto_e + nto_h);

        decompose(tdm.alpha() * 2, nto_e, nto_h);
    }
    else {
    }
 */


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

    if (m_no[1]) {
        out << "NOs (alpha)" << std::endl;
        analysis_p1(out, m_no[0]->get_occ(), nno);
        out << "NOs (beta)" << std::endl;
        analysis_p1(out, m_no[1]->get_occ(), nno);
        analysis_p2(out, m_no[2]->get_occ());
    }
    else {
        vec ev = m_no[0]->get_occ() * 2.0;
        out << "NOs" << std::endl;
        analysis_p1(out, ev, nno);
        analysis_p2(out, ev);
    }
}


void no_analysis::export_orbitals(orbital_printer_i &pr) const {

    if (m_no[1]) {
        orbital_selector s_a, s_b;
        build_selector(m_no[0]->get_occ(), s_a);
        build_selector(m_no[1]->get_occ(), s_b);
        pr.perform(orbital_type::no, *m_no[0], s_a, *m_no[1], s_b);
    }
    else {
        orbital_selector s;
        build_selector(m_no[0]->get_occ(), s);
        pr.perform(orbital_type::no, *m_no[0], s);
    }
}


void no_analysis::analysis_p1(std::ostream &out, const vec &ev, size_t nno) {

    double nelec = accu(ev);
    size_t ihomo = ev.n_elem - (size_t)(nelec + 0.5);

    size_t min = (ihomo > nno ? ihomo - nno : 0);
    size_t max = (ihomo + nno < ev.n_elem ? ihomo + nno : ev.n_elem);

    std::string offset(2, ' ');
    out << offset << "Occupation of frontier NOs:";
    out << std::fixed << std::setprecision(4);
    for (size_t i = min; i < max; i++) out << " " << std::setw(6) << ev(i);
    out << std::endl;
    out << offset << "Number of electrons: ";
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
    out << "Number of unpaired electrons: n_u = ";
    out << std::setw(8) << nu;
    out << ", n_u,nl = ";
    out << std::setw(8) << nunl << std::endl;

    out << "NO participation ratio (PR_NO): ";
    out << std::setw(9) << std::setprecision(6) << (nu * nu) / nu2 << std::endl;
}


void no_analysis::build_selector(const vec &e, orbital_selector &sel) {

    size_t ntot = e.size(), nelec = accu(e);
    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    sel.select(true, ntot - nelec, ntot, 1, true);
    sel.select(false, ntot - 2 * nelec, ntot - nelec, 1, true);
}

} // namespace libwfa


