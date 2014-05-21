#include "pop_mulliken.h"

namespace libwfa {

pop_mulliken::pop_mulliken(const std::vector<size_t> &b2p,
        const arma::Mat<double> &s) : m_nparts(0), m_b2p(b2p), m_s(s) {

    // Compute the number of atoms according to the provided map:
    //      largest atom number in map + 1
    for (std::vector<size_t>::const_iterator i = m_b2p.begin();
            i != m_b2p.end(); i++) {

        m_nparts = std::max(m_nparts, *i);
    }
    m_nparts++;
}


void pop_mulliken::perform(
        const arma::Mat<double> &d_bb, std::vector<double> &p) const {

    p.clear();
    p.resize(m_nparts, 0.0);

    for (size_t i = 0; i != m_b2p.size(); i++) {
        p[m_b2p[i]] += dot(d_bb.row(i), m_s.row(i));
    }
}

} // namespace libwfa


