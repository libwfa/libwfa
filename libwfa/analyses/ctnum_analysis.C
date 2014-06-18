#include "ctnum_analysis.h"

namespace libwfa {

using namespace arma;


ctnum_analysis::ctnum_analysis(const Col<size_t> &b2p) :
    m_nparts(0), m_b2p(b2p) {

    m_nparts = b2p.max() + 1;
}


void ctnum_analysis::perform(const Mat<double> &om_ao, Mat<double> &om) const {

    om.resize(m_nparts, m_nparts);
    om.fill(0.0);

    for (size_t i = 0; i < m_b2p.size(); i++) {

        size_t iat = m_b2p[i];
        for (size_t j = 0; j < m_b2p.size(); j++) {

            size_t jat = m_b2p[j];
            om(iat, jat) += om_ao(i, j);
        }
    }
}


} // namespace libwfa
