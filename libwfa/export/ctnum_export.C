#include <iomanip>
#include <fstream>
#include "ctnum_export.h"

namespace libwfa {

using namespace arma;


void ctnum_export::perform(const ab_matrix &om) const {

    std::string fname(m_prefix + ".om");
    mat om_trace = om.sp_trace();
    do_export(fname, om_trace);

}


void ctnum_export::do_export(const std::string &fname,
    const mat &ct) const {

    std::ofstream out;
    out.open(fname.c_str(), std::ios_base::app);

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
