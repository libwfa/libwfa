#include "pop_mulliken.h"

namespace libwfa {

pop_mulliken::pop_mulliken(const std::vector<size_t> &b2c,
        const arma::Mat<double> &s) : m_natoms(0), m_b2c(b2c), m_s(s) {

    // Compute the number of atoms according to the provided map:
    //      largest atom number in map + 1
    for (std::vector<size_t>::const_iterator i = m_b2c.begin();
            i != m_b2c.end(); i++) {

        m_natoms = std::max(m_natoms, *i);
    }
    m_natoms++;
}


void pop_mulliken::perform(
        const arma::Mat<double> &d_bb, std::vector<double> &p) const {

    p.clear();
    p.resize(m_natoms, 0.0);

    for (size_t i = 0; i != m_b2c.size(); i++) {
        p[m_b2c[i]] += dot(d_bb.row(i), m_s.row(i));
    }
}

} // namespace libwfa


