#include "ctnum_export.h"
#include <iomanip>
#include <fstream>

namespace libwfa {

using namespace arma;

void ctnum_export::perform(const ctnum_data &ct) {
    
    for (ctnum_data::iterator i = ct.begin(); i != ct.end(); i++) {

        std::string fname(m_prefix + ct.name(i) + ".om");

        std::ofstream out;
        out.open(fname.c_str());

        // Set precision and format for doubles
        out << std::setprecision(m_prec) << std::fixed;

        // First header line
        out << ct.name(i) << " ";
        out << ct.energy(i) << " " << ct.osc_strength(i) << std::endl;

        const Mat<double> &data = ct.data(i);
        out << "2 " << data.n_cols << " " << data.n_rows << std::endl;
        // TODO: check which is the correct order for export as a linear array
        for (size_t i1 = 0, j = 0; i1 < data.n_cols; i1++) {
            for (size_t i2 = 0; i2 < data.n_rows; i2++, j++) {
                out << std::setw(m_colwidth) << data(i2, i1);
                if ((j + 1) % m_ncols == 0) out << std::endl;
            }
        }
        out << std::endl;

        out.close();
    }
}

} // namespace libwfa
