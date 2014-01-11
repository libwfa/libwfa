#include "ctnum_analysis.h"

namespace libwfa {

using namespace arma;


ctnum_analysis::ctnum_analysis(const std::vector<size_t> &b2c) :
    m_natoms(0), m_b2c(b2c) {

    // Compute the number of atoms according to the provided map:
    //      largest atom number in map + 1
    for (std::vector<size_t>::const_iterator i = m_b2c.begin();
            i != m_b2c.end(); i++) {

        m_natoms = std::max(m_natoms, *i);
    }
    m_natoms++;
}
        

void ctnum_analysis::perform(const Mat<double> &om_ao, Mat<double> &om_at) {

    om_at.resize(m_natoms, m_natoms);
    om_at.fill(0.0);

    for (size_t i = 0; i < m_b2c.size(); i++) {

        size_t iat = m_b2c[i];
        for (size_t j = 0; j < m_b2c.size(); j++) {

            size_t jat = m_b2c[j];
            om_at(iat, jat) += om_ao(i, j);
        }
    }
}


} // namespace libwfa
