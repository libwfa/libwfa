#include <libwfa/libwfa_exception.h>
#include "density_printer_cube.h"

namespace libwfa {

using namespace arma;


const char density_printer_cube::k_clazz[] = "density_printer_cube";


void density_printer_cube::perform(density_type type, const ab_matrix &dm) {

    if (! m_dt.test(type)) return;

    bool ab_sep = false;
    size_t spin_tr_d = 1;

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    if (dm.is_alpha_eq_beta()) {
        m_core.perform(name, desc, dm.alpha() * 2.0);
    }
    else {
        const mat &dm_a = dm.alpha(),  &dm_b = dm.beta();
        if (ab_sep) {
            m_core.perform(name + "_a", desc + " (alpha part)", dm_a);
            m_core.perform(name + "_b", desc + " (beta part)", dm_b);
        }
        if (spin_tr_d >= 1) {
            mat dm_tr = dm.alpha() + dm.beta();
            m_core.perform(name + "_tot", desc + " (total)", dm_tr);
        }
        if (spin_tr_d >= 2) {
            mat dm_d = dm.alpha() - dm.beta();
            m_core.perform(name + "_sd", desc + " (spin-difference)", dm_d);
        }
    }
}


} // namespace libwfa

