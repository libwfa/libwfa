#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "ev_printer_ndo.h"

namespace libwfa {

using namespace arma;

const char ev_printer_ndo::k_clazz[] = "ev_printer_ndo";


size_t ev_printer_ndo::perform(density_type type,
    const ab_vector &ni, std::ostream &out) const {

    static const char *method =
            "perform(density_type, const ab_vector &, std::ostream &) const";

    if (type != density_type::difference)
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    std::string title("NDOs");
    if (ni.is_alpha_eq_beta()) {
        out << " " << title << std::endl;
        return print(ni.alpha(), out);
    }
    else {
        out << " " << title << "(alpha)" << std::endl;
        size_t na = print(ni.alpha(), out);
        out << " " << title << "(beta)" << std::endl;
        size_t nb = print(ni.beta(), out);
        return std::max(na, nb);
    }
}


size_t ev_printer_ndo::print(const Col<double> &ni, std::ostream &out) const {

    out << std::setw(7) << std::setprecision(4) << std::fixed;
    out << "Leading detachment eigenvalues: ";
    for (size_t i = 0; i < m_nndo; i++)
        out << ni[i];
    out << std::endl;

    out << "Leading attachment eigenvalues: ";
    for (size_t i = 0, j = ni.n_elem - 1; i < m_nndo; i++, j--)
        out << ni[j];
    out << std::endl;

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

    out << "Number of detached / attached electrons: p_D = " << nd;
    out << ", p_A = " << na << std::endl;
    out << std::setw(9) << std::setprecision(6);
    out << "Number of involved orbitals: PR_D = " << nd * nd / nd2;
    out << ", PR_A = " << na * na / na2 << std::endl;

    return ni.n_elem;
}


} // namespace libwfa


