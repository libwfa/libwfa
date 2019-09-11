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
    double ef = type.export_factor();
    bool do_esp = m_dt_esp.test(type);

    if (dm.is_alpha_eq_beta()) {
        m_core.perform(name, desc, dm.alpha() * 2.0 * ef, do_esp);
    }
    else {
        const mat &dm_a = dm.alpha(),  &dm_b = dm.beta();
        if (m_dt_tot.test(type)) {
            m_core.perform(name + "_tot",
                    desc + " (total)", (dm.alpha() + dm.beta()) * ef, do_esp);
            m_core.perform(name + "_sd",
                    desc + " (spin)", dm.alpha() - dm.beta());
        }
        else {
            m_core.perform(name + "_a", desc + " (alpha part)", dm_a * ef, do_esp);
            m_core.perform(name + "_b", desc + " (beta part)", dm_b * ef, do_esp);
        }
    }
}


} // namespace libwfa
