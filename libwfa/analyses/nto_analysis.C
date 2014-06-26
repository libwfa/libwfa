#include <algorithm>
#include <libwfa/core/transformations_dm.h>
#include "nto_analysis.h"

namespace libwfa {

using namespace arma;


void nto_analysis_basic::perform(export_data_i &opr, std::ostream &out) const {

    // Diagonalize particle density matrix
    ab_vector ee;
    ab_matrix ue;
    diagonalize_dm(m_s, m_c, m_edm, ee, ue);
    size_t ne = m_pr.perform(density_type::particle, ee, out);

    // Diagonalize hole density matrix
    ab_vector eh;
    ab_matrix uh;
    diagonalize_dm(m_s, m_c, m_hdm, eh, uh);
    size_t nh = m_pr.perform(density_type::hole, eh, out);

    // Form full matrix u and vector e (properly sorted)
    bool aeqb = ue.is_alpha_eq_beta();
    ab_matrix c_nto(aeqb);
    ab_vector n_nto(aeqb);
    ab_orbital_selector s_nto(aeqb);

    size_t nhtot = uh.alpha().n_cols, ntot = nhtot + ue.alpha().n_cols;
    s_nto.alpha() = orbital_selector(ntot);
    s_nto.alpha().select(true, nhtot - nh, nhtot - 1, 1);
    s_nto.alpha().select(false, nhtot, nhtot + ne, 1);

    c_nto.alpha() = join_rows(uh.alpha(), fliplr(ue.alpha()));
    n_nto.alpha() = join_cols(eh.alpha() * -1., flipud(ee.alpha()));

    if (! aeqb) {
        nhtot = uh.beta().n_cols;
        ntot = nhtot + ue.beta().n_cols;
        s_nto.beta() = orbital_selector(ntot);
        s_nto.beta().select(true, nhtot - nh, nhtot - 1, 1);
        s_nto.beta().select(false, nhtot, nhtot + ne, 1);

        c_nto.beta() = join_rows(uh.beta(), fliplr(ue.beta()));
        n_nto.beta() = join_cols(eh.beta() * -1., flipud(ee.beta()));
    }

    opr.perform(orbital_type::nto, c_nto, n_nto, s_nto);
}


void nto_analysis::perform(ab_matrix &edm, ab_matrix &hdm,
    export_data_i &opr, std::ostream &out) const {

    form_eh(m_s, m_tdm, edm, hdm);

    nto_analysis_basic(m_s, m_c, edm, hdm, m_pr).perform(opr, out);
}


} // namespace libwfa


