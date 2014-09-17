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

    if (ni.is_alpha_eq_beta()) {
        out << "NDOs:" << std::endl;
        size_t n = print(ni.alpha() * 2., out);
        out << std::endl;
        return n;
    }
    else {
        out << "NDOs (alpha):" << std::endl;
        size_t na = print(ni.alpha(), out);
        out << "NDOs (beta):" << std::endl;
        size_t nb = print(ni.beta(), out);
        out << std::endl;
        return std::max(na, nb);
    }
}


size_t ev_printer_ndo::print(const vec &ni, std::ostream &out) const {

    // Compute # attached and detached electrons first 
    double na = 0.0, na2 = 0.0, nd = 0.0, nd2 = 0.0;
    size_t i = 0, nndo = 0;
    for (; i < ni.n_elem && ni(i) < 0; i++, nndo++) {
        double cur = ni(i);
        nd += cur; nd2 += cur * cur;
    }
    nd *= -1;
    for (; i < ni.n_elem; i++) {
        double cur = ni(i);
        na += cur; na2 += cur * cur;
    }
    nndo = std::min(nndo, ni.n_elem - nndo);
    nndo = std::min(nndo, m_nndo);
 
    std::string offset(2, ' ');
    out << std::setprecision(4) << std::fixed;
    out << offset << "Leading detachment eigenvalues: ";
    for (size_t i = 0; i < nndo; i++) 
        out << std::setw(9) << ni(i);
    out << std::endl;

    out << offset << "Leading attachment eigenvalues: ";
    for (size_t i = 0, j = ni.n_elem - 1; i < nndo; i++, j--) 
        out << std::setw(9) << ni(j);
    out << std::endl;

    out << offset << "Number of detached / attached electrons: p_D = ";
    out << std::setw(7) << nd;
    out << ", p_A = " << std::setw(7) << na << std::endl;
    out << std::setprecision(6);
    out << offset << "Number of involved orbitals: PR_D = ";
    out << std::setw(9) << (nd * nd) / nd2;
    out << ", PR_A = ";
    out << std::setw(9) << (na * na) / na2 << std::endl;

    return nndo;
}


} // namespace libwfa


