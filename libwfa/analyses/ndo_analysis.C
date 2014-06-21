#include <libwfa/core/transformations_dm.h>
#include "ndo_analysis.h"

namespace libwfa {

using namespace arma;


void ndo_analysis::perform(ab_matrix &at, ab_matrix &de,
    export_data_i &opr, std::ostream &out) const {

    ab_matrix u;
    ab_vector ev;
    diagonalize_dm(m_c, m_ddm, ev, u);

    size_t nndo = m_pr.perform(density_type::difference, ev, out);

    // Form full matrix u and vector e (properly sorted)

    bool aeqb = u.is_alpha_eq_beta();
    ab_orbital_selector s(aeqb);

    size_t ntot = ev.alpha().n_elem;
    s.alpha() = orbital_selector(ntot);
    s.alpha().select(true, 0, nndo);
    s.alpha().select(false, ntot - nndo, ntot - 1);
    if (! aeqb) {
        ntot = ev.beta().n_elem;
        s.beta() = orbital_selector(ntot);
        s.beta().select(true, 0, nndo);
        s.beta().select(false, ntot - nndo, ntot - 1);
    }

    opr.perform(orbital_type::ndo, u, ev, s);

    // Compute u^-1 = u' * s
    u.alpha() = u.alpha().t() * m_s;
    if (! aeqb) u.beta() = u.beta().t() * m_s;

    form_ad(ev, u, at, de);
}


void ndo_analysis::perform(export_data_i &pr, std::ostream &out) const {

    ab_matrix at, de;
    perform(at, de, pr, out);

    pr.perform(density_type::attach, at);
    pr.perform(density_type::detach, de);
}


} // namespace libwfa
