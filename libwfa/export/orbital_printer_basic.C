#include <libwfa/libwfa_exception.h>
#include "orbital_printer_basic.h"

namespace libwfa {

using namespace arma;


const char orbital_printer_basic::k_clazz[] = "orbital_printer_basic";


void orbital_printer_basic::perform(orbital_type type,
        const orbital_data &orb, const orbital_selector &s) {

    static const char method[] = "perform(orbital_type, "
            "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) return;

    // Do basic error checking
    if (s.n_indexes() != orb.get_coeff().n_cols ||
            s.n_indexes() != orb.get_occ().n_rows) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }

    m_out << m_title << " - " << type << std::endl;

    const vec &ev = orb.get_occ();
    const mat &c = orb.get_coeff();
    Col<uword> el = s.get_selected_arma();

    for (size_t i = 0; i < el.n_elem; i++) {

        std::ostringstream oss;
        oss << " Occupation " << ev(el(i));
        c.col(el(i)).t().print(m_out, oss.str());
    }
}


void orbital_printer_basic::perform(orbital_type type,
        const orbital_data &orb_a, const orbital_selector &s_a,
        const orbital_data &orb_b, const orbital_selector &s_b) {

    static const char method[] = "perform(orbital_type, "
            "const orbital_data &, const orbital_selector &, "
            "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) return;

    // Do basic error checking
    if (s_a.n_indexes() != orb_a.get_coeff().n_cols ||
            s_a.n_indexes() != orb_a.get_occ().n_rows) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }
    if (s_b.n_indexes() != orb_b.get_coeff().n_cols ||
            s_b.n_indexes() != orb_b.get_occ().n_rows) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (beta).");
    }

    m_out << m_title << " - " << type << std::endl;
    m_out << "Alpha spin part:" << std::endl;

    const vec &e_a = orb_a.get_occ();
    const mat &c_a = orb_a.get_coeff();
    Col<uword> el_a = s_a.get_selected_arma();

    for (size_t i = 0; i < el_a.n_elem; i++) {

        std::ostringstream oss;
        oss << " Occupation " << e_a(el_a(i));
        c_a.col(el_a(i)).t().print(m_out, oss.str());
    }

    const vec &e_b = orb_b.get_occ();
    const mat &c_b = orb_b.get_coeff();
    Col<uword> el_b = s_b.get_selected_arma();

    for (size_t i = 0; i < el_b.n_elem; i++) {

        std::ostringstream oss;
        oss << " Occupation " << e_b(el_b(i));
        c_b.col(el_b(i)).t().print(m_out, oss.str());
    }
}


} // namespace libwfa

