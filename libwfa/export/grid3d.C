#include <algorithm>
#include <libwfa/libwfa_exception.h>
#include "grid3d.h"

using namespace arma;

namespace libwfa {


grid3d::grid3d() :
    m_vec(3, 4, fill::zeros), m_npts(3, fill::ones), m_ntotal(1) {

    for(size_t i = 0; i < 3; i++) m_vec(i, i + 1) = 1.0;
}


grid3d::grid3d(unsigned int n, double d) :
    m_vec(3, 4, fill::zeros), m_npts(3, fill::ones), m_ntotal(1) {

    for(size_t i = 0; i < 3; i++) {
        m_npts(i) = n;
        m_ntotal *= n;
        m_vec(i, 0) = - (d * (n - 1)) / 2.;
        m_vec(i, i + 1) = d;
    }
}


grid3d::grid3d(const Col<unsigned int> &n, const vec &d) :
    m_vec(3, 4, fill::zeros), m_npts(n), m_ntotal(1) {

    m_vec.col(0) -= d % (n - 1) / 2.;
    m_vec.cols(1, 3) = diagmat(d);
    for(unsigned int i = 0; i < 3; i++) m_ntotal *= n(i);
}


void grid3d::set_origin(const vec &x0) {

    m_vec.col(0) = x0;
}


void grid3d::set_direction(unsigned int i,
    unsigned int ni, const vec &ei) {
#ifdef LIBWFA_DEBUG
    check_idx(i);
#endif

    m_ntotal = m_ntotal / m_npts(i) * ni;
    m_npts(i) = ni;
    m_vec.col(i + 1) = ei;
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


size_t grid3d::build_pts(size_t i0, mat &pts) const {

    if (pts.n_rows != 3) {
        libwfa_exception("grid3d", "build_pts(size_t, arma::mat &)",
            __FILE__, __LINE__, "pts");
    }

    // Compute index triple for start position and adjust size
    Col<unsigned int> ipos(3, fill::zeros);
    size_t ii = i0, sz = pts.n_cols;
    for (unsigned int i = 0, j = 2; i < 3; i++, j--) {
        ipos(j) = ii % m_npts(j);
        ii /= m_npts(j);
    }
    sz = std::min(sz, m_ntotal - i0);

    if (sz == 0) return sz;

    // Compute first element
    pts.col(0) = m_vec.col(0) + m_vec.cols(1, 3) * ipos;

    if (sz == 1) return sz;

    // Compute increments
    mat inc(3, 3);
    for (unsigned int  i = 0; i < 3; i++) {
       subview<double> ci = inc.col(i);
       ci = m_vec.col(i + 1);
       for (unsigned int j = i + 1; j < 3; j++)
           ci -= m_vec.col(j + 1) * (m_npts(j) - 1);
    }

    // Loop over all subsequent elements
    for (size_t i = 1; i < sz; i++) {
        // Determine increment
        unsigned int j = 2;
        for (; j > 0; j--) {
            ipos(j)++;
            if (ipos(j) < m_npts(j)) break;
            ipos(j) = 0;
        }
        // Add increment to previous element
        pts.col(i) = pts.col(i - 1) + inc.col(j);
    }

    return sz;
}


void grid3d::check_idx(unsigned int idx) {

    if (idx >= 3) {
        throw libwfa_exception("grid3d", "check_idx()",
                    __FILE__, __LINE__, "idx");
     }
}


} // namespace libwfa



