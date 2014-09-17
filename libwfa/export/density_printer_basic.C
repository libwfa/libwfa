#include <libwfa/libwfa_exception.h>
#include "density_printer_basic.h"

namespace libwfa {

using namespace arma;


const char density_printer_basic::k_clazz[] = "density_printer_basic";


void density_printer_basic::perform(density_type type, const ab_matrix &dm) {

    if (! m_dt.test(type)) return;

    m_out << m_title << " - " << type << std::endl;
    if (dm.is_alpha_eq_beta()) {
        dm.alpha().print(m_out);
    }
    else {
        dm.alpha().print(m_out, "Alpha spin part:");
        dm.beta().print(m_out, "Beta spin part:");
    }
}


} // namespace libwfa

