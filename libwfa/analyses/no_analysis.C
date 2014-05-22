#include <algorithm>
#include <libwfa/core/transformations_dm.h>
#include "no_analysis.h"

namespace libwfa {

using namespace arma;


void no_analysis::perform(const ab_matrix &sdm,
        export_orbitals_i &pr_o, std::ostream &out) const {

    ab_matrix c_no;
    ab_vector n_no;
    diagonalize_dm(m_c, sdm, n_no, c_no);

    bool aeqb = c_no.is_alpha_eq_beta();

    // Perform spin-traced calculation
    if (! aeqb) {
        ab_matrix sdm2(true);
        sdm2.alpha() = sdm.alpha() + sdm.beta();

        ab_matrix c2_no;
        ab_vector n2_no;
        diagonalize_dm(m_c, sdm2, n2_no, c2_no);

        n2_no.alpha() *= 0.5;

        m_pr.perform(density_type::state, n2_no, out);

        ab_selector s2_no(true);
        s2_no.alpha().select_all();

        pr_o.perform(orbital_type::no, c2_no, n2_no, s2_no);
    }

    m_pr.perform(density_type::state, n_no, out);

    // Form full matrix u and vector e (properly sorted)


    ab_selector s_no(aeqb);
    s_no.alpha().select_all();
    if (! aeqb) s_no.beta().select_all();

    pr_o.perform(orbital_type::no, c_no, n_no, s_no);
}


} // namespace libwfa


