#include <libwfa/libwfa_exception.h>
#include "export_data_print.h"

namespace libwfa {

using namespace arma;


const char export_data_print::k_clazz[] = "export_data_print";


void export_data_print::perform(density_type type, const ab_matrix &dm) {

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


void export_data_print::perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_orbital_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    if (! m_ot.test(type)) return;

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

    const orbital_selector &sa = s.alpha();
    const Col<double> &eva = ev.alpha();
    const Mat<double> &ca = coeff.alpha();
    Col<uword> ela = sa.get_selected_arma();

    for (size_t i = 0; i < ela.n_elem; i++) {

        std::ostringstream oss;
        oss << " Occupation " << eva(ela(i));
        ca.col(ela(i)).t().print(m_out, oss.str());
    }

    if (aeqb) return;

    const orbital_selector &sb = s.beta();
    const Col<double> &evb = ev.beta();
    const Mat<double> &cb = coeff.beta();
    Col<uword> elb = sb.get_selected_arma();

    for (size_t i = 0; i < elb.n_elem; i++) {

        std::ostringstream oss;
        oss << " Occupation " << evb(elb(i));
        cb.col(elb(i)).t().print(m_out, oss.str());
    }
}


} // namespace libwfa

