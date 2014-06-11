#include <algorithm>
#include <libwfa/core/transformations_dm.h>
#include "nto_analysis.h"

namespace libwfa {

using namespace arma;


void nto_analysis::perform(const ab_matrix_pair &dm, ab_matrix_pair &u,
    export_data_i &opr, std::ostream &out) const {

    // Diagonalize particle density matrix
    ab_vector ee;
    diagonalize_dm(m_c, dm.first, ee, u.first);
    size_t ne = m_pr.perform(density_type::particle, ee, out);

    // Diagonalize hole density matrix
    ab_vector eh;
    diagonalize_dm(m_c, dm.second, eh, u.second);
    size_t nh = m_pr.perform(density_type::hole, eh, out);

    // Form full matrix u and vector e (properly sorted)
    bool aeqb = u.first.is_alpha_eq_beta();
    ab_matrix c_nto(aeqb);
    ab_vector n_nto(aeqb);
    ab_selector s_nto(aeqb);

    size_t ntot_a = u.first.alpha().n_cols + u.second.alpha().n_cols;
    s_nto.alpha() = selector(ntot_a);
    s_nto.alpha().select(0, nh);
    s_nto.alpha().select(ntot_a - ne, ntot_a);

    c_nto.alpha() = join_cols(flipud(u.second.alpha()), u.first.alpha());
    n_nto.alpha() = join_cols(flipud(eh.alpha()) * -1., ee.alpha());

    if (! aeqb) {
        size_t ntot_b = u.first.alpha().n_cols + u.second.beta().n_cols;
        s_nto.beta() = selector(ntot_b);
        s_nto.beta().select(0, nh);
        s_nto.beta().select(ntot_b - ne, ntot_b);

        c_nto.beta() = join_cols(flipud(u.second.beta()), u.first.beta());
        n_nto.beta() = join_cols(flipud(eh.beta()) * -1., ee.beta());
    }

    opr.perform(orbital_type::nto, c_nto, n_nto, s_nto);
}


void nto_analysis::perform(const ab_matrix &tdm, ab_matrix_pair &eh,
    export_data_i &opr, std::ostream &out) const {

    form_eh(m_s, tdm, eh.first, eh.second);

    ab_matrix_pair u;
    nto_analysis::perform(eh, u, opr, out);
}


} // namespace libwfa


