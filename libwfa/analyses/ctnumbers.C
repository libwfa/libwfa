#include <libwfa/core/transformations_dm.h>
#include "ctnumbers.h"

namespace libwfa {


using namespace arma;


void ctnumbers::perform(const ab_matrix &tdm,
        ab_matrix &om, double (&om_tot)[2]) const {

    ab_matrix om_ao;
    form_om(m_s, tdm, om_ao);

    om_tot[0] = accu(om_ao.alpha());
    m_analysis.perform(om_ao.alpha(), om.alpha());

    if (tdm.is_alpha_eq_beta()) {
        om_tot[1] = om_tot[0];
    }
    else {
        om_tot[1] = accu(om_ao.beta());
        m_analysis.perform(om_ao.beta(), om.beta());
    }
}


} // namespace libwfa



