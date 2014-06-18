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

    std::string title;
    if (type == density_type::particle)
        title = "Electron:";
    else if (type == density_type::hole)
        title = "Hole:";
    else
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    if (ni.is_alpha_eq_beta()) {
        out << " " << title << std::endl;
        return print(ni.alpha(), out);
    }
    else {
        out << " a-" << title << std::endl;
        size_t na = print(ni.alpha(), out);

        out << " b-" << title << std::endl;
        size_t nb = print(ni.beta(), out);

        return std::max(na, nb);
    }
}


size_t ev_printer_nto::print(const Col<double> &ni, std::ostream &out) const {

    out << "Leading SVs: ";
    out << std::setw(7) << std::setprecision(4) << std::fixed;
    for (size_t i = 0; i < m_nnto; i++) out << ni[i];
    out << std::endl;

    double total = accu(ni);
    out << std::setw(9) << std::setprecision(6) << std::fixed;
    out << "Sum of SVs: " << total << std::endl;
    out << "NTO participation ratio (PR_NTO): " << total * total / dot(ni, ni);
    out << std::endl;
    Col<uword> x = find(ni < m_thresh, 1, "first");

    return x(0);
}


} // namespace libwfa


