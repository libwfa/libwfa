#include "ndo_analysis.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void ndo_analysis::perform(const ab_matrix &ddm,
        export_densities_i &dm_print, export_orbitals_i &ndo_print,
        ev_data_i &pr) const {

    ab_matrix u;
    ab_vector ev;
    diagonalize_dm(m_c, ddm, ev, u);

    size_t n = pr.perform(dm_type::difference, ev);

    // Form full matrix u and vector e (properly sorted)

    bool aeqb = u.is_alpha_eq_beta();
    ab_selector s(aeqb);
    s.alpha().select_all();
    if (! aeqb) s.beta().select_all();

    ndo_print.perform(u, ev, s);


    // Compute u^-1 = u' * s
    u.alpha() = u.alpha().t() * m_s;
    if (! aeqb) u.beta() = u.beta().t() * m_s;

    ab_matrix att, det;
    form_ad(ev, u, att, det);

    dm_print.perform(dm_type::attach, att);
    dm_print.perform(dm_type::detach, det);
}

} // namespace libwfa


