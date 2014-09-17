#include <libwfa/libwfa_exception.h>
#include "pop_loewdin.h"

namespace libwfa {

using namespace arma;


pop_loewdin::pop_loewdin(const mat &s, const uvec &b2p) :
    m_nparts(b2p.max() + 1), m_b2p(b2p) {

    mat u;
    vec e;
    eig_sym(e, u, s);
    m_sh = u * diagmat(sqrt(e)) * u.t();
}


void pop_loewdin::perform(const mat &d_bb, vec &p) const {

    p = vec(m_nparts, fill::zeros);
    for (size_t i = 0; i != m_b2p.size(); i++) {
        p(m_b2p(i)) -= dot(m_sh.col(i), d_bb * m_sh.col(i));
    }
}


} // namespace libwfa


