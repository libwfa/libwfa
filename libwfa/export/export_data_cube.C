#include <libwfa/libwfa_exception.h>
#include "export_data_cube.h"

namespace libwfa {

using namespace arma;


const char export_data_cube::k_clazz[] = "export_data_cube";


void export_data_cube::perform(density_type type, const ab_matrix &dm, bool ab_sep,
        size_t spin_tr_d) {

    if (! m_dt.test(type)) return;

    std::string name, desc;
    name = m_id + "_" + type.convert();
    {
    std::ostringstream ss; ss << m_desc << " " << type;
    desc = ss.str();
    }

    if (dm.is_alpha_eq_beta()) {
        const Mat<double> &dmx = 2. * dm.alpha();
        m_core.perform(name, desc, dmx);
    }
    else {
        const Mat<double> &dm_a = dm.alpha(),  &dm_b = dm.beta();
        if (ab_sep) {
            m_core.perform(name + "_a", desc + " (alpha part)", dm_a);
            m_core.perform(name + "_b", desc + " (beta part)", dm_b);
        }
        if (spin_tr_d >= 1) {
            Mat<double> dm_tr = dm.alpha() + dm.beta();
            m_core.perform(name + "_sp-tr", desc + " (spin-traced)", dm_tr);
        }
        if (spin_tr_d >= 2) {
            Mat<double> dm_d = dm.alpha() - dm.beta();
            m_core.perform(name + "_sp-diff", desc + " (spin-difference)", dm_d);
        }
    }
}


void export_data_cube::perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ene, const ab_orbital_selector &s) {

    static const char method[] = "perform(const ab_matrix &, "
            "const ab_vector &, const ab_selector &)";

    if (! m_ot.test(type)) return;

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
    const std::string &desc, const Mat<double> &c, const orbital_selector &s) {

    Mat<double> cc = c.cols(s.get_selected_arma());

    m_core.perform(name, desc, s.get_selected(), cc);

}


} // namespace libwfa

