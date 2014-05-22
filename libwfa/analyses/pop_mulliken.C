#include <libwfa/libwfa_exception.h>
#include "pop_mulliken.h"

namespace libwfa {


pop_mulliken::pop_mulliken(const arma::Mat<double> &s,
    const std::vector<size_t> &b2p, const std::vector<double> &p0) :
    m_nparts(p0.size()), m_s(s), m_b2p(b2p), m_p0(p0) {

    size_t nparts = 0;
    // Compute the number of atoms according to the provided map:
    //      largest atom number in map + 1
    for (std::vector<size_t>::const_iterator i = m_b2p.begin();
            i != m_b2p.end(); i++) {

        nparts = std::max(nparts, *i);
    }
    nparts++;

    if (m_nparts == 0) m_nparts = nparts;
    else if (m_nparts != nparts) {
        libwfa_exception("pop_mulliken", "pop_mulliken(..)",
                __FILE__, __LINE__, "nparts");
    }
}


void pop_mulliken::perform(
        const arma::Mat<double> &d_bb, std::vector<double> &p) const {

    if (m_p0.size() != m_nparts) {
        p.clear();
        p.resize(m_nparts, 0.0);
    }
    else {
        p = m_p0;
    }

    for (size_t i = 0; i != m_b2p.size(); i++) {
        p[m_b2p[i]] -= dot(d_bb.row(i), m_s.row(i));
    }
}


} // namespace libwfa


