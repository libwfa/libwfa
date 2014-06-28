#include <libwfa/libwfa_exception.h>
#include "pop_mulliken.h"

namespace libwfa {

using namespace arma; 


pop_mulliken::pop_mulliken(const Mat<double> &s, const Col<size_t> &b2p, 
    const Col<double> &p0) :
    m_nparts(p0.n_elem), m_s(s), m_b2p(b2p), m_p0(p0) {

    size_t nparts = m_b2p.max() + 1;

    if (m_nparts == 0) m_nparts = nparts;
    else if (m_nparts != nparts) {
        libwfa_exception("pop_mulliken", "pop_mulliken(..)",
                __FILE__, __LINE__, "nparts");
    }
}


void pop_mulliken::perform(const Mat<double> &d_bb, Col<double> &p) const {

    if (m_p0.n_elem != m_nparts) {
        p = Col<double>(m_nparts, fill::zeros);
    }
    else {
        p = m_p0;
    }

    for (size_t i = 0; i != m_b2p.size(); i++) {
        p(m_b2p(i)) -= dot(d_bb.row(i), m_s.row(i));
    }
}


} // namespace libwfa


