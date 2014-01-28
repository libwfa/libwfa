#include <algorithm>
#include <iomanip>
#include "libwfa_exception.h"
#include "nto_data_print.h"

namespace libwfa {

using namespace arma;

const char nto_data_print::k_clazz[] = "nto_data_print";


size_t nto_data_print::perform(density_type type, const ab_vector &ni) {

    static const char *method = "perform(density_type, const ab_vector &)";

    std::string title;
    if (type == density_type::particle)
        title = "Electron:";
    else if (type == density_type::hole)
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


} // namespace libwfa


