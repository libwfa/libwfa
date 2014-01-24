#include "export_densities_cube.h"
#include "libwfa_exception.h"

namespace libwfa {

using namespace arma;


void export_densities_cube::perform(density_type type, const ab_matrix &dm) {

    std::string name(m_prefix + "_" + type.convert());
    if (dm.is_alpha_eq_beta()) {
        const Mat<double> &dmx = dm.alpha();
        m_core.perform(name, dmx);
    }
    else {
        const Mat<double> &dm_a = dm.alpha(),  &dm_b = dm.beta();
        m_core.perform(name + "_a", dm_a);
        m_core.perform(name + "_b", dm_b);
    }
}


} // namespace libwfa

