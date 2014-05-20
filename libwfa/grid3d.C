#include <algorithm>
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

    for(unsigned int i = 0; i < 12; i++) m_vec[i] = 0.0;
    for(unsigned int i = 0, j = 3; i < 3; i++, j++) {
        m_npts[i] = n[i];
        m_vec[i] = - (d[i] * n[i]) / 2.;
        m_vec[j] = d[i];
    }
}


void grid3d::set_origin(const double (&x0)[3]) {

    for (unsigned int i = 0; i < 3; i++) m_vec[i] = x0[i];
}


void grid3d::set_direction(unsigned int i, unsigned int ni,
        const double (&ei)[3]) {
#ifdef LIBWFA_DEBUG
    check_idx(i);
#endif

    m_npts[i] = ni;
    for (unsigned int j = 0, k = 3 + 3 * i; j < 3; j++, k++) {
        m_vec[k] = ei[j];
    }
}

void grid3d::check() const {

#ifdef LIBWFA_DEBUG
    for (unsigned int i = 0; i < 3; i++) {
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


void grid3d::build_pts(double *pts, size_t i0, size_t sz) const {

    // Compute index triple for start position and adjust size
    size_t npts = 1, ipos[3];
    for (unsigned int i = 0, j = 2; i < 3; i++, j--) {
        npts *= m_npts[j];
        ipos[j] = i0 % npts;
    }
    sz = std::min(sz, npts - i0);

    if (sz == 0) return;

    // Compute first element
    double *p1 = pts;
    for (unsigned int k = 0; k < 3; k++) {
        *p1 = m_vec[k];
        for (unsigned int i = 0, ik = 3 + k; i < 3; i++, ik += 3) {
            *p1 += ipos[i] * m_vec[ik];
        }
        p1++;
    }

    if (sz == 1) return;

    // Compute increments
    double inc[9];
    { // Scope of pi1 and pv
    double *pi1 = inc;
    const double *pv = m_vec + 3;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int k = 0; k < 3; k++) *pi1++ = *pv++;

        double *pi2 = inc;
        for (unsigned int j = 0; j < i; j++) {
            pv -= 3;
            for (unsigned int k = 0; k < 3; k++) {
                *pi2++ -= *pv++ * (m_npts[i] - 1);
            }
        }
    }
    } // End scope of pi1 and pv

    // Loop over all subsequent elements
    const double *p0 = pts;
    for (size_t i = 1; i < sz; i++) {
        // Determine increment
        unsigned int j = 2;
        for (; j > 0; j--) {
            ipos[j]++;
            if (ipos[j] < m_npts[j]) break;
            ipos[j] = 0;
        }
        // Add increment to previous element
        const double *pi = inc + j * 3;
        for (unsigned int k = 0; k < 3; k++) *p1++ = *p0++ + *pi++;
    }
}


void grid3d::check_idx(unsigned int idx) {

    if (idx >= 3) {
        throw libwfa_exception("grid3d", "check_idx()",
                    __FILE__, __LINE__, "idx");
     }
}


} // namespace libwfa



