#include "ctnumbers.h"
#include "transformations_dm.h"

namespace libwfa {


using namespace arma;


void ctnumbers::perform(const state_info &state, const ab_matrix &tdm) {

    ab_matrix om_ao;
    form_om(m_s, tdm, om_ao);

    ctnum_data data;
    if (om_ao.is_alpha_eq_beta()) {
        Mat<double> &om = data.add(state.convert('_'),
                state.energy, state.osc_strength);

        m_omega[0] = m_omega[1] = accu(om_ao.alpha());
        m_analysis.perform(om_ao.alpha(), om);
    }
    else {
        Mat<double> &om_a = data.add(state.convert('_') + "_a",
                state.energy, state.osc_strength);

        m_omega[0] = accu(om_ao.alpha());
        m_analysis.perform(om_ao.alpha(), om_a);

        Mat<double> &om_b = data.add(state.convert('_') + "_b",
                state.energy, state.osc_strength);

        m_omega[1] = accu(om_ao.beta());
        m_analysis.perform(om_ao.beta(), om_b);
    }

    m_printer.perform(data);
}


} // namespace libwfa



