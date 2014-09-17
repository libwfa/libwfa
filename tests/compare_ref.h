#ifndef LIBWFA_COMPARE_REF_H
#define LIBWFA_COMPARE_REF_H

#include <armadillo>
#include <libtest/libtest.h>
#include <libwfa/core/ab_matrix.h>

namespace libwfa {


class compare_ref {
public:
    /** \brief Compare two Armadillo matrices
        \param test Testname
        \param m Matrix
        \param m_ref Reference matrix
        \param thresh Threshold
     **/
    static void compare(const char *test, arma::mat &m,
        arma::mat &m_ref, double thresh)
        throw(libtest::test_exception);

    /** \brief Compare two ab_matrices
        \param test Testname
        \param m Matrix
        \param m_ref Reference matrix
        \param thresh Threshold
     **/
    static void compare(const char *test, ab_matrix &m,
        ab_matrix &m_ref, double thresh)
        throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_COMPARE_REF_H
