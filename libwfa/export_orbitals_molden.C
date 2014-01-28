#include "export_orbitals_molden.h"
#include "libwfa_exception.h"

namespace libwfa {

using namespace arma;


const char export_orbitals_molden::k_clazz[] = "export_orbitals_molden";


export_orbitals_molden::export_orbitals_molden(molden_file_base &file,
    size_t no_a, size_t nv_a, size_t no_b, size_t nv_b) :
    m_file(file) {

    m_norbs[0] = no_a;
    m_norbs[1] = nv_a;
    m_norbs[2] = no_b;
    m_norbs[3] = nv_b;
}


void export_orbitals_molden::perform(orbital_type type, const ab_matrix &coeff,
    const ab_vector &ene, const ab_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    bool aeqb = coeff.is_alpha_eq_beta();

    // Do basic error checking
    if (s.nidx_a() != coeff.ncols_a() || s.nidx_a() != ene.nrows_a()) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }
    if (coeff.ncols_a() != m_norbs[0] + m_norbs[1]) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Sizes (alpha).");
    }
    if (! aeqb) {
        if (s.nidx_b() != coeff.ncols_b() || s.nidx_b() != ene.nrows_b()) {
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                    "Inconsistent sizes (beta).");
        }
        if (coeff.ncols_a() != m_norbs[0] + m_norbs[1]) {
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                    "Sizes (beta).");
        }
    }

    const selector &sa = s.alpha(), &sb = s.beta();
    if (sa.all_selected() && sb.all_selected()) {

        m_file.perform(coeff, ene, m_norbs[0], m_norbs[2]);
        return;
    }

    ab_matrix c(aeqb);
    ab_vector e(aeqb);

    Col<uword> ia = sa.get_selected_arma();
    size_t no_a = 0, no_b = 0;
    for (; no_a < ia.n_rows && ia(no_a) < m_norbs[0]; no_a++) ;

    c.alpha() = coeff.alpha().cols(ia);
    e.alpha() = ene.alpha().rows(ia);

    if (! aeqb) {
        Col<uword> ib = sb.get_selected_arma();
        size_t no_b = 0;
        for (; no_b < ib.n_rows && ib(no_b) < m_norbs[2]; no_b++) ;

        c.beta() = coeff.beta().cols(ib);
        e.beta() = ene.beta().rows(ib);
    }

    m_file.perform(c, e, no_a, no_b);
}


} // namespace libwfa

