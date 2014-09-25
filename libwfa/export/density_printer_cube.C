#include <libwfa/libwfa_exception.h>
#include "density_printer_cube.h"

namespace libwfa {

using namespace arma;


const char density_printer_cube::k_clazz[] = "density_printer_cube";


void density_printer_cube::perform(density_type type, const ab_matrix &dm) {

    if (! m_dt.test(type)) return;

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    if (dm.is_alpha_eq_beta()) {
        if (type == density_type::hole)
            m_core.perform(name, desc, dm.alpha() * -2.0);
        else
            m_core.perform(name, desc, dm.alpha() * 2.0);
    }
    else {
        const mat &dm_a = dm.alpha(),  &dm_b = dm.beta();
        if (m_dt_tot.test(type)) {
            if (type == density_type::hole) {
                m_core.perform(name + "_tot",
                        desc + " (total)", (dm.alpha() + dm.beta()) * -1.);
                m_core.perform(name + "_sd",
                        desc + " (spin)", dm.alpha() - dm.beta());
            }
            else {
                m_core.perform(name + "_tot",
                        desc + " (total)", dm.alpha() + dm.beta());
                m_core.perform(name + "_sd",
                        desc + " (spin)", dm.alpha() - dm.beta());
            }
        }
        else {
            if (type == density_type::hole) {
                m_core.perform(name + "_a", desc + " (alpha part)", dm_a * -1.);
                m_core.perform(name + "_b", desc + " (beta part)", dm_b * -1.);
            }
            else {
                m_core.perform(name + "_a", desc + " (alpha part)", dm_a);
                m_core.perform(name + "_b", desc + " (beta part)", dm_b);
            }
        }
    }
}


} // namespace libwfa

