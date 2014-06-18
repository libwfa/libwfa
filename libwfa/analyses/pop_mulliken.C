#include <libwfa/libwfa_exception.h>
#include "pop_mulliken.h"

namespace libwfa {


pop_mulliken::pop_mulliken(const arma::Mat<double> &s,
    const arma::Col<size_t> &b2p, const arma::Col<double> &p0) :
    m_nparts(p0.n_elem), m_s(s), m_b2p(b2p), m_p0(p0) {

    size_t nparts = m_b2p.max() + 1;

    if (m_nparts == 0) m_nparts = nparts;
    else if (m_nparts != nparts) {
        libwfa_exception("pop_mulliken", "pop_mulliken(..)",
                __FILE__, __LINE__, "nparts");
    }
}


void pop_mulliken::perform(
        const arma::Mat<double> &d_bb, arma::Col<double> &p) const {

    if (m_p0.n_elem != m_nparts) {
        p.clear();
        p.resize(m_nparts, 0.0);
    }
    else {
        p = m_p0;
    }

    for (size_t i = 0; i != m_b2p.size(); i++) {
        p(m_b2p(i)) -= dot(d_bb.row(i), m_s.row(i));
    }
}


} // namespace libwfa


