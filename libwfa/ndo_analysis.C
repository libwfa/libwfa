#include "ndo_analysis.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void ndo_analysis::perform(const ab_matrix &ddm, ab_matrix_pair &ad,
        export_orbitals_i &ndo_print, ev_data_i &pr) const {

    ab_matrix u;
    ab_vector ev;
    diagonalize_dm(m_c, ddm, ev, u);

    size_t n = pr.perform(density_type::difference, ev);

    // Form full matrix u and vector e (properly sorted)

    bool aeqb = u.is_alpha_eq_beta();
    ab_selector s(aeqb);
    s.alpha().select_all();
    if (! aeqb) s.beta().select_all();

    ndo_print.perform(orbital_type::ndo, u, ev, s);

    // Compute u^-1 = u' * s
    u.alpha() = u.alpha().t() * m_s;
    if (! aeqb) u.beta() = u.beta().t() * m_s;

    form_ad(ev, u, ad.first, ad.second);
}



void ndo_analysis::perform(const ab_matrix &ddm,
        export_densities_i &dm_print, export_orbitals_i &ndo_print,
        ev_data_i &pr) const {

    ab_matrix_pair ad;
    perform(ddm, ad, ndo_print, pr);

    dm_print.perform(density_type::attach, ad.first);
    dm_print.perform(density_type::detach, ad.second);
}

} // namespace libwfa


