#include "pop_analysis_ad.h"

namespace libwfa {

using namespace arma;


void pop_analysis_ad::perform(const ab_matrix &attach, const ab_matrix &detach,
    pop_data &pop) const {

    if (attach.is_alpha_eq_beta()) {

        std::vector<double> &h = pop.add("h+");
        m_analysis.perform(detach.alpha(), h);
        std::vector<double> &e = pop.add("e-");
        m_analysis.perform(attach.alpha(), e);
        std::vector<double> &dq = pop.add("Del q");
        dq.resize(m_analysis.size(), 0.0);

        std::vector<double>::const_iterator ih = h.begin(), ie = e.begin();
        std::vector<double>::iterator iq = dq.begin();
        for (; iq != dq.end(); iq++, ie++, ih++) { (*iq) = (*ih) + (*ie); }
    }
    else {
        std::vector<double> &ha = pop.add("h+ (alpha)");
        m_analysis.perform(detach.alpha(), ha);
        std::vector<double> &hb = pop.add("h+ (beta)");
        m_analysis.perform(detach.beta(), hb);
        std::vector<double> &ea = pop.add("e- (alpha)");
        m_analysis.perform(attach.alpha(), ea);
        std::vector<double> &eb = pop.add("e- (beta)");
        m_analysis.perform(attach.alpha(), eb);
    }
}


} // namespace libwfa




