#include <libwfa/libwfa_exception.h>
#include "export_orbitals_molden.h"

namespace libwfa {

using namespace arma;


const char export_orbitals_molden::k_clazz[] = "export_orbitals_molden";


export_orbitals_molden::export_orbitals_molden(export_molden_i &core,
    const std::string &id, const ot_flag &ot) :
    m_core(core), m_id(id), m_ot(ot) {

}


void export_orbitals_molden::perform(orbital_type type, const ab_matrix &coeff,
    const ab_vector &ene, const ab_orbital_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    if (! m_ot.test(type)) return;

    bool aeqb = coeff.is_alpha_eq_beta();

    // Do basic error checking
    if (s.nidx_a() != coeff.ncols_a() || s.nidx_a() != ene.nrows_a()) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }
    if (! aeqb) {
        if (s.nidx_b() != coeff.ncols_b() || s.nidx_b() != ene.nrows_b()) {
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                    "Inconsistent sizes (beta).");
        }
    }

    std::string name(m_id + "_" + type.convert());

    const orbital_selector &sa = s.alpha(), &sb = s.beta();
    size_t no_a = sa.n_selected(true);
    size_t no_b = sb.n_selected(true);

    ab_matrix c(aeqb);
    ab_vector e(aeqb);

    Col<uword> ia = sa.get_selected_arma();

    c.alpha() = coeff.alpha().cols(ia);
    e.alpha() = ene.alpha().rows(ia);

    if (! aeqb) {
        Col<uword> ib = sb.get_selected_arma();

        c.beta() = coeff.beta().cols(ib);
        e.beta() = ene.beta().rows(ib);
    }

    m_core.perform(name, c, e, no_a, no_b);
}


} // namespace libwfa

