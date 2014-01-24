#include <algorithm>
#include <iomanip>
#include "libwfa_exception.h"
#include "nto_data.h"

namespace libwfa {

using namespace arma;

const char nto_data_print::k_clazz[] = "nto_data_print";


size_t nto_data_print::perform(dm_type type, const ab_vector &ni) {

    static const char *method = "perform(dm_type, const ab_vector &)";

    std::string title;
    if (type == dm_type::particle)
        title = "Electron:";
    else if (type == dm_type::hole)
        title = "Hole:";
    else
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    if (ni.is_alpha_eq_beta()) {
        m_out << " " << title << std::endl;
        return print(ni.alpha());
    }
    else {
        m_out << " a-" << title << std::endl;
        size_t na = print(ni.alpha());

        m_out << " b-" << title << std::endl;
        size_t nb = print(ni.beta());

        return std::max(na, nb);
    }
}


size_t nto_data_print::print(const Col<double> &ni) {

    m_out << "Leading SVs: ";
    for (size_t i = 0; i < m_nnto; i++) {
        m_out << std::setw(7) << std::setprecision(4) << std::fixed << ni[i];
    }
    double total = accu(ni);
    m_out << "Sum of SVs: ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed << total;
    m_out << "NTO participation ratio (PR_NTO): ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed;
    m_out << total * total / dot(ni, ni);
    Col<uword> x = find(ni < m_thresh, 1, "first");

    return x(0);
}


size_t nto_data_extract::perform(dm_type type, const ab_vector &ni) {

    static const char *method = "perform(dm_type, const ab_vector &)";

    if (type != dm_type::particle && type != dm_type::hole) {
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");
    }

    m_sets.push_back(ab_ntoinfo(ni.is_alpha_eq_beta()));

    const Col<double> &ni_a = ni.alpha();
    ntoinfo r_a = m_sets.back().alpha();

    r_a.type = type;
    r_a.total = accu(ni_a);
    r_a.pr = r_a.total * r_a.total / dot(ni_a, ni_a);
    r_a.ni.resize(m_nnto, 0.0);
    for (size_t i = 0; i < m_nnto; i++) r_a.ni[i] = ni_a[i];

    Col<uword> x_a = find(ni_a < m_thresh, 1, "first");

    if (ni.is_alpha_eq_beta()) return x_a(0);

    const Col<double> &ni_b = ni.beta();
    ntoinfo r_b = m_sets.back().beta();

    r_b.type = type;
    r_b.total = accu(ni_b);
    r_b.pr = r_b.total * r_b.total / dot(ni_b, ni_b);
    r_b.ni.resize(m_nnto, 0.0);
    for (size_t i = 0; i < m_nnto; i++) r_b.ni[i] = ni_b[i];

    Col<uword> x_b = find(ni_b < m_thresh, 1, "first");

    return std::max(x_a(0), x_b(0));
}


} // namespace libwfa


