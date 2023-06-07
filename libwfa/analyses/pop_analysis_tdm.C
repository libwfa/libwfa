#include "pop_analysis_tdm.h"

namespace libwfa {

using namespace arma;


void pop_analysis_tdm::perform(pop_data &pop) const {

    vec &ch = pop.add("Trans. (e)");
    m_analysis.perform(m_tdm.alpha() + m_tdm.beta(), ch);

}


} // namespace libwfa
