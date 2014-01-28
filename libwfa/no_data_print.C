#include <algorithm>
#include <iomanip>
#include "libwfa_exception.h"
#include "no_data_print.h"

namespace libwfa {

using namespace arma;

const char no_data_print::k_clazz[] = "no_data_print";


size_t no_data_print::perform(density_type type, const ab_vector &ni) {

    static const char *method = "perform(density_type, const ab_vector &)";

    if (type != density_type::state)
        throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "type.");

    std::string title("NOs");
    if (ni.is_alpha_eq_beta()) {
        m_out << " " << title << "(spin-traced)" << std::endl;
        return print_total(ni.alpha() * 2.0);
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

    double nelec = accu(ni);
    size_t ihomo = ni.n_elem - (nelec + 0.5);

    m_out << "Occupation of frontier NOs: ";
    for (size_t i = 0, j = ihomo - m_nno; i < 2 * m_nno; i++, j++) {
        m_out << std::setw(7) << std::setprecision(4) << std::fixed << ni[j];
    }
    m_out << std::endl;
    m_out << "Number of electrons: ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed << nelec;
    m_out << std::endl;

    return nelec + 0.5;
}

size_t no_data_print::print_total(const Col<double> &ni) {

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

    m_out << "Occupation of frontier NOs: ";
    for (size_t i = 0, j = ihomo - m_nno; i < 2 * m_nno; i++, j++) {
        m_out << std::setw(7) << std::setprecision(4) << std::fixed << ni[j];
    }
    m_out << std::endl;
    m_out << "Number of electrons: ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed << nelec;

    m_out << "Number of unpaired electrons: n_u = ";
    m_out << std::setw(8) << std::setprecision(5) << std::fixed << nu;
    m_out << ", n_u,nl = ";
    m_out << std::setw(8) << std::setprecision(5) << std::fixed << nunl;
    m_out << std::endl;
    m_out << "NO participation ratio (PR_NO): ";
    m_out << std::setw(9) << std::setprecision(6) << std::fixed;
    m_out << nu * ni / nu2 << std::endl;

    return nelec / 2. + 0.5;
}


} // namespace libwfa


