#include <libwfa/core/transformations_dm.h>
#include "ndo_analysis.h"

namespace libwfa {

using namespace arma;


void ndo_analysis::perform(ab_matrix &at, ab_matrix &de,
    ab_matrix &u, ab_vector &ev,
    export_data_i &opr, std::ostream &out) const {

    diagonalize_dm(m_s, m_c, m_ddm, ev, u);

    size_t nndo = m_pr.perform(density_type::difference, ev, out);

    // Form full matrix u and vector e (properly sorted)

    bool aeqb = u.is_alpha_eq_beta();
    ab_orbital_selector s(aeqb);

    size_t ntot = ev.alpha().n_elem;
    nndo = std::min(2 * nndo, ntot) / 2;
    s.alpha() = orbital_selector(ntot);
    s.alpha().select(true, 0, nndo, 1, true);
    s.alpha().select(false, ntot - nndo, ntot, 1, true);
    if (! aeqb) {
        ntot = ev.beta().n_elem;
        s.beta() = orbital_selector(ntot);
        s.beta().select(true, 0, nndo, 1, true);
        s.beta().select(false, ntot - nndo, ntot, 1, true);
    }

    opr.perform(orbital_type::ndo, u, ev, s);

    form_ad(ev, u, at, de);
}

void ndo_analysis::perform(ab_matrix &at, ab_matrix &de,
    export_data_i &opr, std::ostream &out) const {
        
    ab_matrix u;
    ab_vector ev;
    perform(at, de, u, ev, opr, out);        
}

void ndo_analysis::perform(export_data_i &opr, std::ostream &out) const {

    ab_matrix at, de;
    ab_matrix u;
    ab_vector ev;
    perform(at, de, u, ev, opr, out);

    opr.perform(density_type::attach, at);
    opr.perform(density_type::detach, de);
}


} // namespace libwfa
