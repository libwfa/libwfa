#include <libwfa/libwfa_exception.h>
#include "export_densities_print.h"

namespace libwfa {

using namespace arma;


void export_densities_print::perform(density_type type, const ab_matrix &dm) {

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

