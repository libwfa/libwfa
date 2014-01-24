#include <algorithm>
#include <iomanip>
#include "libwfa_exception.h"
#include "no_data.h"

namespace libwfa {

using namespace arma;

const char no_data_print::k_clazz[] = "no_data_print";


size_t no_data_print::perform(dm_type type, const ab_vector &ni) {

    static const char *method = "perform(dm_type, const ab_vector &)";

    if (type != dm_type::state)
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    std::string title("NOs");
    if (ni.is_alpha_eq_beta()) {
        m_out << " " << title << "(spin-traced)" << std::endl;
        return print(ni.alpha() * 2.0);
    }
    else {
        m_out << " " << title << "(alpha)" << std::endl;
        size_t na = print(ni.alpha());

        m_out << " " << title << "(beta)" << std::endl;
        size_t nb = print(ni.beta());

        return std::max(na, nb);
    }
}


size_t no_data_print::print(const Col<double> &ni) {

    m_out << "Leading detachment eigenvalues: ";
    for (size_t i = 0; i < m_nno; i++) {
        m_out << std::setw(7) << std::setprecision(4) << std::fixed << ni[i];
    }
    m_out << "Leading attachment eigenvalues: ";
    for (size_t i = 0, j = ni.n_elem - 1; i < m_nno; i++, j--) {
        m_out << std::setw(7) << std::setprecision(4) << std::fixed << ni[j];
    }

    double na = 0.0, na2 = 0.0, nd = 0.0, nd2 = 0.0;

    size_t i = 0;
    for (; i < ni.n_elem && ni(i) < 0; i++) {
        double cur = ni(i);
        nd += cur; nd2 += cur * cur;
    }
    nd *= -1;
    for (; i < ni.n_elem; i++) {
        double cur = ni(i);
        na += cur; na2 += cur * cur;
    }

    m_out << "Number of detached / attached electrons: p_D = ";
    m_out << std::setw(7) << std::setprecision(4) << std::fixed << nd;
    m_out << ", p_A = ";
    m_out << std::setw(7) << std::setprecision(4) << std::fixed << na;
    m_out << "Number of involved orbitals: PR_D = ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed;
    m_out << nd * nd / nd2;
    m_out << ", PR_A = ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed;
    m_out << na * na / na2;

    return ni.n_elem;
}


} // namespace libwfa


