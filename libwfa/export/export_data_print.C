#include <libwfa/libwfa_exception.h>
#include "export_data_print.h"

namespace libwfa {

using namespace arma;


const char export_data_print::k_clazz[] = "export_data_print";


void export_data_print::perform(density_type type, const ab_matrix &dm) {

    m_out << m_title << " - " << type << std::endl;
    if (dm.is_alpha_eq_beta()) {
        dm.alpha().print(m_out);
    }
    else {
        dm.alpha().print(m_out, "Alpha spin part:");
        dm.beta().print(m_out, "Beta spin part:");
    }
}


void export_data_print::perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    bool aeqb = coeff.is_alpha_eq_beta();

    // Do basic error checking
    if (s.nidx_a() != coeff.ncols_a() || s.nidx_a() != ev.nrows_a()) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                "Inconsistent sizes (alpha).");
    }
    if (! aeqb) {
        if (s.nidx_b() != coeff.ncols_b() || s.nidx_b() != ev.nrows_b()) {
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__,
                    "Inconsistent sizes (beta).");
        }
    }

    m_out << m_title << " - " << type << std::endl;
    if (! aeqb) m_out << "Alpha spin part:";

    const selector &sa = s.alpha();
    const arma::Col<double> &eva = ev.alpha();
    const arma::Mat<double> &ca = coeff.alpha();

    for (size_t i = 0; i < sa.n_indexes(); i++) {
        if (! sa.is_selected(i)) continue;

        std::ostringstream oss;
        oss << " Occupation " << eva[i];
        ca.col(i).t().print(m_out, oss.str());
    }

    if (aeqb) return;

    const selector &sb = s.beta();
    const arma::Col<double> &evb = ev.beta();
    const arma::Mat<double> &cb = coeff.beta();

    for (size_t i = 0; i < sb.n_indexes(); i++) {
        if (! sb.is_selected(i)) continue;

        std::ostringstream oss;
        oss << " Occupation " << evb[i];
        cb.col(i).t().print(m_out, oss.str());
    }
}


} // namespace libwfa

