#include <algorithm>
#include "no_analysis.h"

namespace libwfa {

using namespace arma;


void no_analysis::perform(const ev_printer_i &evpr, export_data_i &opr,
                          std::ostream &out) const {

    bool aeqb = m_u.is_alpha_eq_beta();

    if (! aeqb)
        evpr.perform(density_type::state, m_est, out);    
    
    size_t nelec = evpr.perform(density_type::state, m_e, out);

    // Form full matrix u and vector e (properly sorted)

    ab_orbital_selector s_no(aeqb);

    size_t ntot = m_e.alpha().size();
    s_no.alpha() = orbital_selector(ntot);
    s_no.alpha().select(true, ntot - nelec, ntot, 1, true);
    s_no.alpha().select(false, ntot - 2 * nelec , ntot - nelec, 1, true);
    if (! aeqb) {
        ntot = m_e.beta().size();
        s_no.beta() = orbital_selector(ntot);
        s_no.beta().select(true, ntot - nelec, ntot, 1, true);
        s_no.beta().select(false, ntot - 2 * nelec , ntot - nelec, 1, true);
    }

    opr.perform(orbital_type::no, m_u, m_e, s_no);
}


} // namespace libwfa


