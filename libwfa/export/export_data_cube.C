#include <libwfa/libwfa_exception.h>
#include "export_data_cube.h"

namespace libwfa {

using namespace arma;


const char export_data_cube::k_clazz[] = "export_data_cube";


void export_data_cube::perform(density_type type, const ab_matrix &dm) {

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    if (dm.is_alpha_eq_beta()) {
        const Mat<double> &dmx = dm.alpha();
        m_core.add(name, desc, dmx);
    }
    else {
        const Mat<double> &dm_a = dm.alpha(),  &dm_b = dm.beta();
        m_core.add(name + "_a", desc + " (alpha part)", dm_a);
        m_core.add(name + "_b", desc + " (beta part)", dm_b);
    }
}


void export_data_cube::perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ene, const ab_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    if (s.nidx_a() != coeff.ncols_a() || s.nidx_b() != coeff.ncols_b()) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "sizes");
    }

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    if (coeff.is_alpha_eq_beta()) {
        perform(name, desc, coeff.alpha(), s.alpha());
    }
    else {
        perform(name + "_a", desc + " (alpha part)", coeff.alpha(), s.alpha());
        perform(name + "_b", desc + " (beta part)", coeff.beta(), s.beta());
    }

}


void export_data_cube::perform(const std::string &name,
    const std::string &desc, const Mat<double> &c, const selector &s) {

    if (s.all_selected()) {

        m_core.add(name, desc, s.get_selected(), c);
    }
    else if (! s.none_selected()) {
        Mat<double> cc = c.cols(s.get_selected_arma());

        m_core.add(name, desc, s.get_selected(), cc);
    }

}


} // namespace libwfa

