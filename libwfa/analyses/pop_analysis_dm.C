#include "pop_analysis_dm.h"

namespace libwfa {

using namespace arma;


void pop_analysis_dm::perform(pop_data &pop) const {

    std::vector<double> &ch = pop.add("Charge (e)");
    m_analysis.perform(m_sdm.alpha() + m_sdm.beta(), ch);

    if (! m_sdm.is_alpha_eq_beta()) {

        std::vector<double> &sp = pop.add("Spin (e)");
        m_analysis.perform(m_sdm.alpha() - m_sdm.beta(), sp);
    }
}


} // namespace libwfa




