#include "pop_analysis_dm.h"

namespace libwfa {

using namespace arma;


void pop_analysis_dm::perform(const ab_matrix &sdm, pop_data &pop) const {

    std::vector<double> &ch = pop.add("Charge (e)");
    m_analysis.perform(sdm.alpha() + sdm.beta(), ch);

    if (sdm.is_alpha_eq_beta()) {

        std::vector<double> &sp = pop.add("Spin (e)");
        m_analysis.perform(sdm.alpha() - sdm.beta(), sp);
    }
}


} // namespace libwfa




