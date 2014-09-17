#include "ctnumbers.h"

namespace libwfa {


using namespace arma;


void ctnumbers::perform(ab_matrix &om, double (&om_tot)[2]) const {

    ab_matrix om_ao;
    form_om(m_s, m_tdm, om_ao);

    om_tot[0] = accu(om_ao.alpha());
    m_analysis.perform(om_ao.alpha(), om.alpha());

    if (m_tdm.is_alpha_eq_beta()) {
        om.set_alpha_eq_beta();
        om_tot[1] = om_tot[0];
    }
    else {
        om.set_alpha_neq_beta();
        om_tot[1] = accu(om_ao.beta());
        m_analysis.perform(om_ao.beta(), om.beta());
    }
}


void ctnumbers::form_om(const mat &s, const ab_matrix &tdm, ab_matrix &om) {

    if (tdm.is_alpha_eq_beta()) {
        om.set_alpha_eq_beta();
        form_om(s, tdm.alpha(), om.alpha());
    }
    else {
        om.set_alpha_neq_beta();
        form_om(s, tdm.alpha(), om.alpha());
        form_om(s, tdm.beta(), om.beta());
    }
}


} // namespace libwfa



