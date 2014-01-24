#include <algorithm>
#include "nto_data.h"

namespace libwfa {

using namespace arma;


size_t nto_data_print::perform(dm_type type, const ab_vector &ni) {

    const Col<double> &ni_a = ni.alpha();

    double tot_a = accu(ni_a);
    double pr_a = tot_a * tot_a / dot(ni_a, ni_a);
    std::vector<double> nx_a(m_nnto, 0.0);
    for (size_t i = 0; i < m_nnto; i++) nx_a[i] = ni_a[i];
    Col<uword> x_a = find(ni_a < m_thresh, 1, "first");

    if (ni.is_alpha_eq_beta()) return x_a(0);

    const Col<double> &ni_b = ni.beta();

    double tot_b = accu(ni_b);
    double pr_b = tot_b * tot_b / dot(ni_b, ni_b);
    std::vector<double> nx_b(m_nnto, 0.0);
    for (size_t i = 0; i < m_nnto; i++) nx_b[i] = ni_b[i];
    Col<uword> x_b = find(ni_b < m_thresh, 1, "first");

    return std::max(x_a(0), x_b(0));
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


