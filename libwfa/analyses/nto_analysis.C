#include <algorithm>
#include "nto_analysis.h"

namespace libwfa {

using namespace arma;


void nto_analysis_basic::perform(const ev_printer_i &evpr,
        export_data_i &opr, std::ostream &out) const {

    size_t ne = evpr.perform(density_type::particle, m_ee, out);
    size_t nh = evpr.perform(density_type::hole, m_eh, out);

    // Form full matrix u and vector e (properly sorted)
    bool aeqb = m_ue.is_alpha_eq_beta();
    ab_matrix c_nto(aeqb);
    ab_vector n_nto(aeqb);
    ab_orbital_selector s_nto(aeqb);

    size_t nhtot = m_uh.alpha().n_cols, ntot = nhtot + m_ue.alpha().n_cols;
    s_nto.alpha() = orbital_selector(ntot);
    s_nto.alpha().select(true, nhtot - nh, nhtot, 1);
    s_nto.alpha().select(false, nhtot, nhtot + ne, 1);

    c_nto.alpha() = join_rows(m_uh.alpha(), fliplr(m_ue.alpha()));
    n_nto.alpha() = join_cols(m_eh.alpha() * -1., flipud(m_ee.alpha()));

    if (! aeqb) {
        nhtot = m_uh.beta().n_cols;
        ntot = nhtot + m_ue.beta().n_cols;
        s_nto.beta() = orbital_selector(ntot);
        s_nto.beta().select(true, nhtot - nh, nhtot, 1);
        s_nto.beta().select(false, nhtot, nhtot + ne, 1);
        c_nto.beta() = join_rows(m_uh.beta(), fliplr(m_ue.beta()));
        n_nto.beta() = join_cols(m_eh.beta() * -1., flipud(m_ee.beta()));
    }

    opr.perform(orbital_type::nto, c_nto, n_nto, s_nto);
}

void nto_analysis::perform(ab_matrix &edm, ab_matrix &hdm,
    export_data_i &opr, std::ostream &out) const {

    form_eh(m_s, m_tdm, edm, hdm);

    nto_analysis_basic(m_s, m_c, edm, hdm).perform(m_pr, opr, out);
}


} // namespace libwfa


