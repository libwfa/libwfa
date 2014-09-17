#include <iomanip>
#include <fstream>
#include "ctnum_export.h"

namespace libwfa {

using namespace arma;


void ctnum_export::perform(const ab_matrix &om,
    const double (&om_tot)[2], std::ostream &out) const {
    
    size_t w = 10, prec = 4;
    out << "  omega = ";
    out << std::setprecision(prec) << std::fixed;
    out << std::setw(w) << om_tot[0] + om_tot[1] << " (alpha: ";
    out << std::setw(w) << om_tot[0] << ", beta: ";
    out << std::setw(w) << om_tot[1] << ")" << std::endl << std::endl;

    if (om.is_alpha_eq_beta()) {
        std::string fname(m_prefix + ".om");
        do_export(fname, om.alpha());
    }
    else {
        std::string fname_a(m_prefix + "_a.om");
        do_export(fname_a, om.alpha());

        std::string fname_b(m_prefix + "_b.om");
        do_export(fname_b, om.beta());
    }
}


void ctnum_export::do_export(const std::string &fname,
    const mat &ct) const {

    std::ofstream out;
    out.open(fname.c_str());

    // Set precision and format for doubles
    out << std::setprecision(m_prec) << std::fixed;

    // First header line
    out << m_desc << std::endl;

    out << "2 " << ct.n_cols << " " << ct.n_rows << std::endl;
    // TODO: check which is the correct order for export as a linear array
    for (size_t i1 = 0, j = 0; i1 < ct.n_cols; i1++) {
        for (size_t i2 = 0; i2 < ct.n_rows; i2++, j++) {
            out  << std::setw(m_colwidth) << ct(i2, i1);
            if ((j + 1) % m_ncols == 0) out << std::endl;
        }
    }
    out << std::endl;

    out.close();
}

} // namespace libwfa
