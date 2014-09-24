#include "ctnum_analysis.h"

namespace libwfa {

using namespace arma;


ctnum_analysis::ctnum_analysis(const mat &s, const uvec &b2p) :
    m_nparts(0), m_s(s), m_b2p(b2p) {

    m_nparts = b2p.max() + 1;
}


void ctnum_analysis::perform(const mat &tdm, mat &om) const {

    mat om_ao;
    form_om(m_s, tdm, om_ao);

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


void ctnum_analysis::form_om(const arma::mat &s,
    const arma::mat &tdm, arma::mat &om) {

    om = 0.5 * ((tdm * s) % (s * tdm) + tdm % (s * tdm * s));
}


} // namespace libwfa
