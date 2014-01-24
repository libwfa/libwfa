#include <algorithm>
#include "nto_analysis.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void nto_analysis::perform(const ab_matrix_pair &dm, ab_matrix_pair &u,
    export_orbitals_i &nto_print, nto_data_i &pr) const {

    ab_vector eh, ee;
    diagonalize_dm(m_c, dm.first, ee, u.first);
    diagonalize_dm(m_c, dm.second, eh, u.second);

    size_t ne = pr.perform(dm_type::edm, ee);
    size_t nh = pr.perform(dm_type::hdm, eh);

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

    nto_print.perform(c_nto, n_nto, s_nto);
}


void nto_analysis::perform(const ab_matrix &tdm, ab_matrix_pair &av,
    export_densities_i &dm_print, export_orbitals_i &nto_print,
    nto_data_i &pr) const {

    ab_matrix_pair dm;
    form_eh(m_s, tdm, dm.first, dm.second);

    dm_print.perform(dm_type::edm, dm.first);
    dm_print.perform(dm_type::hdm, dm.second);

    ab_matrix_pair u;
    nto_analysis::perform(dm, u, nto_print, pr);

    av.first.alpha() += dm.first.alpha();
    av.second.alpha() += dm.second.alpha();

    if (tdm.is_alpha_eq_beta()) return;

    av.first.beta()  += dm.first.beta();
    av.second.beta() += dm.second.beta();
}


} // namespace libwfa


