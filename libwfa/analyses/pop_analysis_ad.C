#include "pop_analysis_ad.h"

namespace libwfa {

using namespace arma;


void pop_analysis_ad::perform(pop_data &pop) const {

    if (m_at.is_alpha_eq_beta()) {

        Col<double> &h = pop.add("h+");
        m_analysis.perform(m_de.alpha() * -2., h);
        Col<double> &e = pop.add("e-");
        m_analysis.perform(m_at.alpha() *  2., e);
        Col<double> &dq = pop.add("Del q");
        dq = h + e;
    }
    else {
        Col<double> &ha = pop.add("h+ (alpha)");
        m_analysis.perform(m_de.alpha() * -1., ha);
        Col<double> &hb = pop.add("h+ (beta)");
        m_analysis.perform(m_de.beta() * -1., hb);
        Col<double> &ea = pop.add("e- (alpha)");
        m_analysis.perform(m_at.alpha(), ea);
        Col<double> &eb = pop.add("e- (beta)");
        m_analysis.perform(m_at.beta(), eb);
    }
}


} // namespace libwfa




