#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "ev_printer_nto.h"

namespace libwfa {

using namespace arma;

const char ev_printer_nto::k_clazz[] = "ev_printer_nto";


size_t ev_printer_nto::perform(density_type type,
    const ab_vector &ni, std::ostream &out) const {

    static const char *method =
            "perform(density_type, const ab_vector &, std::ostream &)";

    bool aeqb = ni.is_alpha_eq_beta();
    std::string pt;
    if (type == density_type::particle)
        pt = "electron";
    else if (type == density_type::hole)
        pt = "hole";
    else
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    if (ni.is_alpha_eq_beta()) {
        out << "NTOs (" << pt << "):" << std::endl;
        size_t n = print(ni.alpha(), out);
        out << std::endl;
        return n;
    }
    else {
        out << "NTOs (" << pt << ", alpha):" << std::endl;
        size_t na = print(ni.alpha(), out);
        out << "NTOs (" << pt << ", beta):" << std::endl;
        size_t nb = print(ni.beta(), out);
        out << std::endl;

        return std::max(na, nb);
    }
}


size_t ev_printer_nto::print(const Col<double> &ni, std::ostream &out) const {

    std::string offset(2, ' ');
    out << offset << "Leading SVs: ";
    out << std::setprecision(4) << std::fixed;
    for (size_t i = 0, j = ni.n_rows - 1; i < m_nnto; i++, j--) 
        out << std::setw(9) << ni(j);
    out << std::endl;

    double total = accu(ni);
    out << std::setprecision(6) << std::fixed;
    out << offset << "Sum of SVs:  " << std::setw(9) << total << std::endl;
    out << offset << "Participation ratio (PR_NTO): " 
        << std::setw(9) << total * total / dot(ni, ni);
    out << std::endl;
    Col<uword> x = find(ni > m_thresh, 1, "first");

    return x(0);
}


} // namespace libwfa


