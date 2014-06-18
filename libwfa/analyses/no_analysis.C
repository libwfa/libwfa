#include <algorithm>
#include <libwfa/core/transformations_dm.h>
#include "no_analysis.h"

namespace libwfa {

using namespace arma;


void no_analysis::perform(export_data_i &opr, std::ostream &out) const {

    ab_matrix c_no;
    ab_vector n_no;
    diagonalize_dm(m_c, m_sdm, n_no, c_no);

    bool aeqb = c_no.is_alpha_eq_beta();

    // Perform spin-traced calculation
    if (! aeqb) {
        ab_matrix sdm2(true);
        sdm2.alpha() = m_sdm.alpha() + m_sdm.beta();

        ab_matrix c2_no;
        ab_vector n2_no;
        diagonalize_dm(m_c, sdm2, n2_no, c2_no);

        n2_no.alpha() *= 0.5;

        m_pr.perform(density_type::state, n2_no, out);
    }

    size_t nelec = m_pr.perform(density_type::state, n_no, out);

    // Form full matrix u and vector e (properly sorted)

    ab_orbital_selector s_no(aeqb);

    size_t ntot_a = n_no.alpha().size();
    s_no.alpha() = orbital_selector(ntot_a);
    s_no.alpha().select(true, ntot_a - nelec, ntot_a - 1, 1, true);
    s_no.alpha().select(false, ntot_a - 2 * nelec , ntot_a - nelec - 1, 1, true);
    if (! aeqb) {
        size_t ntot_b = n_no.beta().size();
        s_no.beta() = orbital_selector(ntot_b);
        s_no.beta().select(true, ntot_b - nelec, ntot_b - 1, 1, true);
        s_no.beta().select(false, ntot_b - 2 * nelec , ntot_b - nelec - 1, 1, true);
    }

    opr.perform(orbital_type::no, c_no, n_no, s_no);
}


} // namespace libwfa


