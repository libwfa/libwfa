#include <libwfa/core/transformations_dm.h>
#include "ndo_analysis.h"

namespace libwfa {

using namespace arma;


void ndo_analysis::perform(const ab_matrix &ddm, ab_matrix_pair &ad,
        export_orbitals_i &pr_o, ev_printer_i &pr_e) const {

    ab_matrix u;
    ab_vector ev;
    diagonalize_dm(m_c, ddm, ev, u);

    pr_e.perform(density_type::difference, ev);

    // Form full matrix u and vector e (properly sorted)

    bool aeqb = u.is_alpha_eq_beta();
    ab_selector s(aeqb);
    s.alpha().select_all();
    if (! aeqb) s.beta().select_all();

    pr_o.perform(orbital_type::ndo, u, ev, s);

    // Compute u^-1 = u' * s
    u.alpha() = u.alpha().t() * m_s;
    if (! aeqb) u.beta() = u.beta().t() * m_s;

    form_ad(ev, u, ad.first, ad.second);
}



void ndo_analysis::perform(const ab_matrix &ddm,
        export_densities_i &pr_d, export_orbitals_i &pr_o,
        ev_printer_i &pr_e) const {

    ab_matrix_pair ad;
    perform(ddm, ad, pr_o, pr_e);

    pr_d.perform(density_type::attach, ad.first);
    pr_d.perform(density_type::detach, ad.second);
}

} // namespace libwfa


