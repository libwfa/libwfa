#include "pop_analysis_ad.h"

namespace libwfa {

using namespace arma;


void pop_analysis_ad::perform(const ab_matrix_pair &ad, pop_print_i &pr) const {

    pop_data res;

    if (ad.first.is_alpha_eq_beta()) {

        std::vector<double> &h = res.add("h+");
        m_analysis.perform(ad.second.alpha(), h);
        std::vector<double> &e = res.add("e-");
        m_analysis.perform(ad.first.alpha(), e);
        std::vector<double> &dq = res.add("Del q");
        dq.resize(m_analysis.size(), 0.0);

        std::vector<double>::const_iterator ih = h.begin(), ie = e.begin();
        std::vector<double>::iterator iq = dq.begin();
        for (; iq != dq.end(); iq++, ie++, ih++) { (*iq) = (*ih) + (*ie); }
    }
    else {
        std::vector<double> &ha = res.add("h+ (alpha)");
        m_analysis.perform(ad.second.alpha(), ha);
        std::vector<double> &hb = res.add("h+ (beta)");
        m_analysis.perform(ad.second.beta(), hb);
        std::vector<double> &ea = res.add("e- (alpha)");
        m_analysis.perform(ad.first.alpha(), ea);
        std::vector<double> &eb = res.add("e- (beta)");
        m_analysis.perform(ad.first.alpha(), eb);
    }
    pr.perform(res);
}


} // namespace libwfa




