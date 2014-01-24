#include <algorithm>
#include <iomanip>
#include "nto_data.h"

namespace libwfa {

using namespace arma;


size_t nto_data_print::perform(dm_type type, const ab_vector &ni) {

    if (ni.is_alpha_eq_beta()) {
        if (type == dm_type::hdm) m_out << " Hole:";
        else if (type == dm_type::edm) { m_out << " Electron:"; }
        m_out << std::endl;

        return print(ni.alpha());
    }
    else {
        if (type == dm_type::hdm) m_out << " a-Hole:";
        else if (type == dm_type::edm) m_out << " a-Electron:";
        m_out << std::endl;
        size_t na = print(ni.alpha());

        if (type == dm_type::hdm) m_out << " b-Hole:";
        else if (type == dm_type::edm) m_out << " b-Electron:";
        m_out << std::endl;
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


