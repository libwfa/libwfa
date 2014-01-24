#include "export_densities_cube.h"
#include "libwfa_exception.h"

namespace libwfa {

using namespace arma;


void export_densities_cube::perform(const state_info &si,
        dm_type type, const ab_matrix &dm) {

    std::string base = type.convert() + "." + si.convert('.');

    if (dm.is_alpha_eq_beta()) {
        const Mat<double> &dmx = dm.alpha();
        m_core.perform(base, dmx);
    }
    else {
        const Mat<double> &dm_a = dm.alpha(),  &dm_b = dm.beta();
        m_core.perform(base + ".alpha", dm_a);
        m_core.perform(base + ".beta", dm_b);
    }
}


} // namespace libwfa

