#include <libwfa/libwfa_exception.h>
#include "test_orbital_printer.h"

namespace libwfa {

using namespace arma;


const char test_orbital_printer::k_clazz[] = "test_orbital_printer";


void test_orbital_printer::perform(orbital_type type,
        const orbital_data &orb, const orbital_selector &s) {

    static const char method[] = "perform(orbital_type, "
            "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Orbital type not active.");
    }

    // Do basic error checking
    if (s.n_indexes() != orb.get_coeff().n_cols ||
            s.n_indexes() != orb.get_occ().n_rows) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }

    check(orb.get_occ(), s);
}


void test_orbital_printer::perform(orbital_type type,
        const orbital_data &orb_a, const orbital_selector &s_a,
        const orbital_data &orb_b, const orbital_selector &s_b) {

    static const char method[] = "perform(orbital_type, "
            "const orbital_data &, const orbital_selector &, "
            "const orbital_data &, const orbital_selector &)";

    if (! m_ot.test(type)) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Orbital type not active.");
    }

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

    check(orb_a.get_occ(), s_a);
    check(orb_b.get_occ(), s_b);
}


void test_orbital_printer::check(const vec &e, const orbital_selector &s) {

    static const char method[] =
            "check(const vec &, const orbital_selector &)";

    for (size_t i = 0; i < e.size(); i++) {

        if (s.is_selected(i)) {
            if (fabs(e(i)) < m_thresh)
                throw libwfa_exception(k_clazz, method,
                        __FILE__, __LINE__, "el(i) < thresh");
        }
        else {
            if (fabs(e(i)) > m_thresh)
                throw libwfa_exception(k_clazz, method,
                        __FILE__, __LINE__, "el(i) > thresh");
        }
    }
}


} // namespace libwfa

