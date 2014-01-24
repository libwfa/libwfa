#include "export_orbitals_cube.h"
#include "libwfa_exception.h"

namespace libwfa {

using namespace arma;


void export_orbitals_cube::perform(const ab_matrix &coeff,
        const ab_vector &ene, const ab_selector &s) {

    static const char *clazz = "export_orbitals_cube";
    static const char *method = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    if (s.nidx_a() != coeff.ncols_a() || s.nidx_b() != coeff.ncols_b()) {
        throw libwfa_exception(clazz, method, __FILE__, __LINE__, "sizes");
    }

    const selector &sa = s.alpha();
    const Mat<double> &ca = coeff.alpha();
    if (sa.all_selected()) {

        m_core.perform(cube::data_type::orb_a, sa.get_selected(), ca);
    }
    else if (! sa.none_selected()) {
        Mat<double> ca_copy = ca.cols(sa.get_selected_arma());

        m_core.perform(cube::data_type::orb_a, sa.get_selected(), ca_copy);
    }

    if (coeff.is_alpha_eq_beta()) return;

    const selector &sb = s.beta();
    const Mat<double> &cb = coeff.beta();
    if (sb.all_selected()) {

        m_core.perform(cube::data_type::orb_b, sb.get_selected(), cb);
    }
    else if (! sb.none_selected()) {
        Mat<double> cb_copy = cb.cols(sb.get_selected_arma());

        m_core.perform(cube::data_type::orb_b, sb.get_selected(), cb_copy);
    }
}


} // namespace libwfa

