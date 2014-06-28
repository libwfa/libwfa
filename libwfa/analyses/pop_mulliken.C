#include <libwfa/libwfa_exception.h>
#include "pop_mulliken.h"

namespace libwfa {

using namespace arma;


pop_mulliken::pop_mulliken(const Mat<double> &s, const Col<size_t> &b2p) :
    m_nparts(0), m_s(s), m_b2p(b2p) {

    m_nparts = m_b2p.max() + 1;
}


void pop_mulliken::perform(const Mat<double> &d_bb, Col<double> &p) const {

	p = Col<double>(m_nparts, fill::zeros);
    for (size_t i = 0; i != m_b2p.size(); i++) {
        p(m_b2p(i)) -= dot(d_bb.row(i), m_s.row(i));
    }
}


} // namespace libwfa


