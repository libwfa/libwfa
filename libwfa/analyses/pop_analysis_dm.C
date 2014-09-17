#include "pop_analysis_dm.h"

namespace libwfa {

using namespace arma;


void pop_analysis_dm::perform(pop_data &pop) const {

    vec &ch = pop.add("Charge (e)");
    m_analysis.perform(m_sdm.alpha() + m_sdm.beta(), ch);
    ch += m_p0;

    if (! m_sdm.is_alpha_eq_beta()) {

        vec &sp = pop.add("Spin (e)");
        m_analysis.perform(m_sdm.alpha() - m_sdm.beta(), sp);
    }
}


} // namespace libwfa




