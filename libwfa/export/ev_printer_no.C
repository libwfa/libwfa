#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "ev_printer_no.h"

namespace libwfa {

using namespace arma;

const char ev_printer_no::k_clazz[] = "ev_printer_no";


size_t ev_printer_no::perform(density_type type,
    const ab_vector &ni, std::ostream &out) const {

    static const char *method =
            "perform(density_type, const ab_vector &, ostream &out) const";

    if (type != density_type::state)
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    if (ni.is_alpha_eq_beta()) {
        out << "NOs (spin-traced)" << std::endl;
        size_t n = print_total(ni.alpha() * 2.0, out);
        out << std::endl;
        return n;
    }
    else {
        out << "NOs (alpha)" << std::endl;
        size_t na = print(ni.alpha(), out);
        out << "NOs (beta)" << std::endl;
        size_t nb = print(ni.beta(), out);
        out << std::endl;

        return std::max(na, nb);
    }
}


size_t ev_printer_no::print(const Col<double> &ni, std::ostream &out) const {

    double nelec = accu(ni);
    size_t ihomo = ni.n_elem - (size_t)(nelec + 0.5);

    std::string offset(2, ' ');
    out << offset << "Occupation of frontier NOs: ";
    out << std::fixed << std::setprecision(4);
    size_t min = (ihomo > m_nno ? ihomo - m_nno : 0);
    size_t max = (ihomo + m_nno <= ni.n_elem ? ihomo + m_nno :  ni.n_elem);
    for (size_t i = min; i < max; i++) out << std::setw(9) << ni(i);
    out << std::endl;
    out << offset << "Number of electrons: "; 
    out << std::setprecision(6) << std::setw(9) << nelec << std::endl;

    return nelec + 0.5;
}

size_t ev_printer_no::print_total(const Col<double> &ni, 
    std::ostream &out) const {

    double nelec = 0.0, nu = 0.0, nu2 = 0.0, nunl = 0.0;
    for (size_t i = 0; i < ni.n_elem; i++) {
        double n = ni(i), nn = 2. - n;
        double tmp = std::min(n, nn);
        nelec += n;
        nu += tmp;
        nu2 += tmp * tmp;
        nunl += n * n * nn * nn;
    }
    size_t ihomo = ni.n_elem - (size_t)(nelec + 0.5) / 2;
    size_t imin = (ihomo > m_nno) ? ihomo - m_nno : 0;
    size_t imax = (ihomo + m_nno <= ni.n_elem) ? ihomo + m_nno : ni.n_elem;
    std::string offset(2, ' ');
    out << std::setprecision(4) << std::fixed;
    out << offset << "Occupation of frontier NOs: ";
    for (size_t i = ihomo, j = imax - 1; i < imax; i++, j--) 
        out << std::setw(9) << ni(j);
    out << " | ";
    for (size_t i = imin, j = ihomo - 1; i < ihomo; i++, j--) 
        out << std::setw(9) << ni(j);
    out << std::endl;

    out << offset << "Number of electrons: ";
    out << std::setprecision(6) << std::setw(9) << nelec << std::endl;

    out << std::setprecision(5);
    out << offset << "Number of unpaired electrons: n_u = ";
    out << std::setw(8) << nu;
    out << ", n_u,nl = ";
    out << std::setw(8) << nunl << std::endl;

    out << offset << "NO participation ratio (PR_NO): ";
    out << std::setw(9) << std::setprecision(6) << (nu * nu) / nu2 << std::endl;

    return nelec / 2. + 0.5;
}


} // namespace libwfa


