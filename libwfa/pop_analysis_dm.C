#include "pop_analysis_dm.h"

namespace libwfa {

using namespace arma;


void pop_analysis_dm::perform(const ab_matrix &sdm, pop_print_i &pr) const {

    pop_data res;
    std::vector<double> &ch = res.add("Charge (e)");
    m_analysis.perform(sdm.alpha() + sdm.beta(), ch);

    if (sdm.is_alpha_eq_beta()) {

        std::vector<double> &sp = res.add("Spin (e)");
        m_analysis.perform(sdm.alpha() - sdm.beta(), sp);
    }
    pr.perform(res);
}


} // namespace libwfa




