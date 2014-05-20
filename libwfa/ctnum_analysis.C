#include "ctnum_analysis.h"

namespace libwfa {

using namespace arma;


ctnum_analysis::ctnum_analysis(const std::vector<size_t> &b2p) :
    m_nparts(0), m_b2p(b2p) {

    // Compute the number of parts according to the provided map:
    //      largest part number in map + 1
    for (std::vector<size_t>::const_iterator i = m_b2p.begin();
            i != m_b2p.end(); i++) {

        m_nparts = std::max(m_nparts, *i);
    }
    m_nparts++;
}


void ctnum_analysis::perform(const Mat<double> &om_ao,
        Mat<double> &om_at) const {

    om_at.resize(m_nparts, m_nparts);
    om_at.fill(0.0);

    for (size_t i = 0; i < m_b2p.size(); i++) {

        size_t iat = m_b2p[i];
        for (size_t j = 0; j < m_b2p.size(); j++) {

            size_t jat = m_b2p[j];
            om_at(iat, jat) += om_ao(i, j);
        }
    }
}


} // namespace libwfa
