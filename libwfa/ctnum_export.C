#include "ctnum_export.h"
#include <iomanip>
#include <fstream>

namespace libwfa {

using namespace arma;


void ctnum_export::perform(const ab_matrix &ct) {
    
    if (ct.is_alpha_eq_beta()) {

        std::string fname(m_prefix + ".om");
        do_export(fname, ct.alpha());
    }
    else {

        std::string fname_a(m_prefix + "_a.om");
        do_export(fname_a, ct.alpha());

        std::string fname_b(m_prefix + "_b.om");
        do_export(fname_b, ct.beta());
    }
}


void ctnum_export::do_export(const std::string &fname, const Mat<double> &ct) {

    std::ofstream out;
    out.open(fname.c_str());

    // Set precision and format for doubles
    out << std::setprecision(m_prec) << std::fixed;

    // First header line
    out << m_prefix << " ";
    out << m_energy << " " << m_osc_strength << std::endl;

    out << "2 " << ct.n_cols << " " << ct.n_rows << std::endl;
    // TODO: check which is the correct order for export as a linear array
    for (size_t i1 = 0, j = 0; i1 < ct.n_cols; i1++) {
        for (size_t i2 = 0; i2 < ct.n_rows; i2++, j++) {
            out << std::setw(m_colwidth) << ct(i2, i1);
            if ((j + 1) % m_ncols == 0) out << std::endl;
        }
    }
    out << std::endl;

    out.close();
}

} // namespace libwfa
