#include <libwfa/libwfa_exception.h>
#include "orbital_printer_molden.h"

namespace libwfa {

using namespace arma;


const char orbital_printer_molden::k_clazz[] = "orbital_printer_molden";


orbital_printer_molden::orbital_printer_molden(export_molden_i &core,
    const std::string &id, const ot_flag &ot) :
    m_core(core), m_id(id), m_ot(ot) {

}


void orbital_printer_molden::perform(orbital_type type,
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

    std::string name(m_id + "_" + type.convert());

    size_t no = s.n_selected(true);

    uvec idx = s.get_selected_arma();
    if (idx.n_rows != 0) {
        mat c = orb.get_coeff().cols(idx);
        vec e = orb.get_occ().rows(idx);
        m_core.perform(name, c, e, no);
    }

}


void orbital_printer_molden::perform(orbital_type type,
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

    std::string name(m_id + "_" + type.convert());

    size_t no_a = s_a.n_selected(true);
    size_t no_b = s_b.n_selected(true);

    uvec idx;
    idx = s_a.get_selected_arma();
    mat c_a, c_b;
    mat e_a, e_b;
    if (idx.n_rows != 0) {
        c_a = orb_a.get_coeff().cols(idx);
        e_a = orb_a.get_occ().rows(idx);
    }
    idx = s_b.get_selected_arma();
    if (idx.n_rows != 0) {
        c_b = orb_b.get_coeff().cols(idx);
        e_b = orb_b.get_occ().rows(idx);
    }

    m_core.perform(name, c_a, e_a, no_a, c_b, e_b, no_b);
}


} // namespace libwfa

