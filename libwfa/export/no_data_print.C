#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "no_data_print.h"

namespace libwfa {

using namespace arma;

const char no_data_print::k_clazz[] = "no_data_print";


size_t no_data_print::perform(density_type type,
    const ab_vector &ni, std::ostream &out) const {

    static const char *method =
            "perform(density_type, const ab_vector &, ostream &out) const";

    if (type != density_type::state)
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    std::string title("NOs");
    if (ni.is_alpha_eq_beta()) {
        out << " " << title << " (spin-traced)" << std::endl;
        return print_total(ni.alpha() * 2.0, out);
    }
    else {
        out << " " << title << " (alpha)" << std::endl;
        size_t na = print(ni.alpha(), out);

        out << " " << title << " (beta)" << std::endl;
        size_t nb = print(ni.beta(), out);

        return std::max(na, nb);
    }
}


size_t no_data_print::print(const Col<double> &ni, std::ostream &out) const {

    double nelec = accu(ni);
    size_t ihomo = ni.n_elem - (nelec + 0.5);

    out << "Occupation of frontier NOs: ";
    out << std::fixed << std::setw(7) << std::setprecision(4);
    for (size_t i = 0, j = ihomo - m_nno; i < 2 * m_nno; i++, j++)
        out << ni[j];
    out << std::endl;
    out << std::setw(9) << std::setprecision(6);
    out << "Number of electrons: " << nelec << std::endl;

    return nelec + 0.5;
}

size_t no_data_print::print_total(
    const Col<double> &ni, std::ostream &out) const {

    double nelec = 0.0, nu = 0.0, nu2 = 0.0, nunl = 0.0;
    for (size_t i = 0; i < ni.n_elem; i++) {
        double n = ni(i), nn = 2. - n;
        double tmp = std::min(n, nn);
        nelec += n;
        nu += tmp;
        nu2 += tmp * tmp;
        nunl += n * n * nn * nn;
    }
    size_t ihomo = ni.n_elem - (nelec + 0.5);

    out << std::setw(7) << std::setprecision(4) << std::fixed;
    out << "Occupation of frontier NOs: ";
    for (size_t i = 0, j = ihomo - m_nno; i < 2 * m_nno; i++, j++)
        out << ni[j];
    out << std::endl;

    out << std::setw(9) << std::setprecision(6);
    out << "Number of electrons: " << nelec;

    out << std::setw(8) << std::setprecision(5);
    out << "Number of unpaired electrons: n_u = " << nu;
    out << ", n_u,nl = " << nunl << std::endl;

    out << std::setw(9) << std::setprecision(6);
    out << "NO participation ratio (PR_NO): " << nu * ni / nu2 << std::endl;

    return nelec / 2. + 0.5;
}


} // namespace libwfa


