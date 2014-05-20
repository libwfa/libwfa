#include "grid3d.h"
#include "libwfa_exception.h"

namespace libwfa {


grid3d::grid3d() {
    for(size_t i = 0; i < 12; i++) m_vec[i] = 0.0;
    for(size_t i = 0, j = 3; i < 3; i++, j += 4) {
        m_npts[i] = 1;
        m_vec[j] = 1.0;
    }
}


grid3d::grid3d(unsigned int n, double d) {
    for(size_t i = 0; i < 12; i++) m_vec[i] = 0.0;
    for(size_t i = 0, j = 3; i < 3; i++, j += 4) {
        m_npts[i] = n;
        m_vec[i] = - (d * n) / 2.;
        m_vec[j] = d;
    }
}


grid3d::grid3d(const unsigned int (&n)[3], const double (&d)[3]) {

    for(size_t i = 0; i < 12; i++) m_vec[i] = 0.0;
    for(size_t i = 0, j = 3; i < 3; i++, j++) {
        m_npts[i] = n[i];
        m_vec[i] = - (d[i] * n[i]) / 2.;
        m_vec[j] = d[i];
    }
}


void grid3d::check() const {

#ifdef LIBWFA_DEBUG
    for (size_t i = 0; i < 3; i++) {
        if (m_npts[i] == 0) {
            throw libwfa_exception("grid3d", "check()",
                    __FILE__, __LINE__, "npts");
        }
    }
    double det = m_vec[0] * (m_vec[4] * m_vec[8] - m_vec[5] * m_vec[7])
            - m_vec[1] * (m_vec[3] * m_vec[8] - m_vec[5] * m_vec[6])
            + m_vec[2] * (m_vec[3] * m_vec[7] - m_vec[4] * m_vec[6]);
    if (det == 0.0) {
        throw libwfa_exception("grid3d", "check()",
                __FILE__, __LINE__, "directions");
    }
#endif
}


void grid3d::check_idx(unsigned int idx) {

    if (idx >= 3) {
        throw libwfa_exception("grid3d", "check_idx()",
                    __FILE__, __LINE__, "idx");
     }
}


} // namespace libwfa



