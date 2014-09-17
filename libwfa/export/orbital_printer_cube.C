#include <libwfa/libwfa_exception.h>
#include "orbital_printer_cube.h"

namespace libwfa {

using namespace arma;


const char orbital_printer_cube::k_clazz[] = "orbital_printer_cube";


void orbital_printer_cube::perform(orbital_type type,
        const orbital_data &orb, const orbital_selector &s) {

    static const char method[] = "perform(orbital_type, "
        "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) return;

    if (s.n_indexes() != orb.get_coeff().n_cols) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "sizes");
    }

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    perform(name, desc, orb.get_coeff(), s);
}


void orbital_printer_cube::perform(orbital_type type,
        const orbital_data &orb_a, const orbital_selector &s_a,
        const orbital_data &orb_b, const orbital_selector &s_b) {

    static const char method[] = "perform(orbital_type, "
        "const orbital_data &, const orbital_selector &, "
        "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) return;

    if (s_a.n_indexes() != orb_a.get_coeff().n_cols ||
            s_b.n_indexes() != orb_b.get_coeff().n_cols) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "sizes");
    }

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    perform(name + "_a", desc + " (alpha part)", orb_a.get_coeff(), s_a);
    perform(name + "_b", desc + " (beta part)", orb_b.get_coeff(), s_b);
}


void orbital_printer_cube::perform(const std::string &name,
    const std::string &desc, const mat &c, const orbital_selector &s) {

    mat cc = c.cols(s.get_selected_arma());
    m_core.perform(name, desc, s.get_selected(), cc);

}


} // namespace libwfa

