#include "export_densities_cube.h"
#include "libwfa_exception.h"

namespace libwfa {

using namespace arma;


void export_densities_cube::perform(dm_type type, size_t idx,
        const ab_matrix &dm) {

    // Determine increment of index for specific density types
    size_t inc = 0;
    if (type == dm_type::hdm) inc = 100;
    else if (type == dm_type::edm) inc = 200;

    std::vector<size_t> ii(1, idx + inc);
    Mat<double> data = dm.alpha();
    data.reshape(dm.alpha().n_elem, 1);

    m_core.perform(determine_data_type(type, true), ii, data);

    if (dm.is_alpha_eq_beta()) return;

    data = dm.beta();
    data.reshape(dm.beta().n_elem, 1);
    m_core.perform(determine_data_type(type, false), ii, data);
}


void export_densities_cube::perform(const dm_list &lst) {

    static const char *clazz = "export_densities_cube";
    static const char *method = "perform(const dm_list &)";

    if (lst.size() == 0) return;

    // Determine increment of index for specific density types
    size_t inc = 0;
    if (lst.type() == dm_type::hdm) inc = 100;
    else if (lst.type() == dm_type::edm) inc = 200;

    // Determine number of alpha / beta densities to export and check for
    // equal sizes
    size_t na = 0, nb = 0, dima, dimb;
    for (dm_list::iterator i = lst.begin(); i != lst.end(); i++) {

        const ab_matrix &dm = lst.get_density(i);
        if (na == 0) dima = dm.alpha().n_elem;
        else if (dima != dm.alpha().n_elem) {
            throw libwfa_exception(clazz, method, __FILE__, __LINE__, "dima");
        }
        na++;

        if (dm.is_alpha_eq_beta()) continue;

        if (nb == 0) dimb = dm.beta().n_elem;
        else if (dimb != dm.beta().n_elem) {
            throw libwfa_exception(clazz, method, __FILE__, __LINE__, "dimb");
        }
        nb++;
    }

    std::vector<size_t> ilst(na, 0);
    Mat<double> data(dima, na);

    size_t pos = 0;
    for (dm_list::iterator i = lst.begin(); i != lst.end(); i++) {

        const ab_matrix &dm = lst.get_density(i);

        ilst[pos] = lst.get_index(i) + inc;
        data.col(pos) = Mat<double>(dm.alpha().mem, dima, 1);
        pos++;
    }

    m_core.perform(determine_data_type(lst.type(), true), ilst, data);

    if (nb == 0) return;

    ilst.resize(nb, 0);
    data.resize(dimb, nb);

    pos = 0;
    for (dm_list::iterator i = lst.begin(); i != lst.end(); i++) {

        const ab_matrix &dm = lst.get_density(i);
        if (dm.is_alpha_eq_beta()) continue;

        ilst[pos] = lst.get_index(i) + inc;
        data.col(pos) = Mat<double>(dm.beta().mem, dimb, 1);
        pos++;
    }

    m_core.perform(determine_data_type(lst.type(), false), ilst, data);
}


cube::data_type export_densities_cube::determine_data_type(dm_type type,
        bool alpha) {

    if (type == dm_type::sdm) {
        return (alpha ? cube::data_type::sdm_a : cube::data_type::sdm_b);
    }
    else if (type == dm_type::tdm ||
            type == dm_type::edm || type == dm_type::hdm) {
        return (alpha ? cube::data_type::tdm_a : cube::data_type::tdm_b);
    }
    else if (type == dm_type::adm) {
        return (alpha ? cube::data_type::adm_a : cube::data_type::adm_b);
    }
    else if (type == dm_type::ddm) {
        return (alpha ? cube::data_type::ddm_a : cube::data_type::ddm_b);
    }
    else {
        throw libwfa_exception("export_densities_cube", "determine_data_type",
                __FILE__, __LINE__, "Unknown dm_type.");
    }
}


} // namespace libwfa

