#include "pop_analysis_ad.h"

namespace libwfa {

using namespace arma;


void pop_analysis_ad::perform(pop_data &pop) const {

    if (m_at.is_alpha_eq_beta()) {

        vec &h = pop.add("h+");
        m_analysis.perform(m_de.alpha() * 2., h);
        vec &e = pop.add("e-");
        m_analysis.perform(m_at.alpha() *  2., e);
        vec &dq = pop.add("Del q");
        dq = h + e;
    }
    else {
        vec &ha = pop.add("h+ (alpha)");
        m_analysis.perform(m_de.alpha(), ha);
        vec &hb = pop.add("h+ (beta)");
        m_analysis.perform(m_de.beta(), hb);
        vec &ea = pop.add("e- (alpha)");
        m_analysis.perform(m_at.alpha(), ea);
        vec &eb = pop.add("e- (beta)");
        m_analysis.perform(m_at.beta(), eb);
    }
}


} // namespace libwfa




