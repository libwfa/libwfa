#ifndef LIBWFA_AB_MATRIX_H
#define LIBWFA_AB_MATRIX_H

#include <armadillo>
#include "ab_object.h"

namespace libwfa {

/** \brief Container for matrices with spin property

    The container stores two matrices with alpha- and beta-spin and a flag to
    indicate, if both matrices have to be identical (i.e. spin-restricted
    calculations). The case of alpha == beta will be enforced by the container
    by using the same matrix for alpha-spin and beta-spin.

    Mathematically, an ab_matrix behaves like a block diagonal matrix of the
      form:
    \f[
    \mathbf{A} = \left(\begin{array}{cc}
                  \mathbf{A}_{\alpha\alpha} & 0 \\
                          0 & \mathbf{A}_{\beta\beta}
                 \end{array}\right)
    \f]

    \ingroup libwfa
 **/
class ab_matrix : public ab_object<arma::mat> {
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta matrix
     **/
    ab_matrix(bool aeqb = false) : ab_object<arma::mat>(aeqb) { }

    /** \brief Constructor for matrix with alpha == beta
        \param nrows Number of rows
        \param ncols Number of columns

        If ncols is not given (== 0), a symmetric matrix is constructed
     **/
    ab_matrix(size_t nrows, size_t ncols = 0) : ab_object<arma::mat>(true) {

        alpha() = arma::mat(nrows, ncols == 0 ? nrows : ncols);
    }

    /** \brief Constructor for matrix with alpha != beta
        \param nrows_a Number of rows
        \param ncols_a Number of columns
        \param nrows_b Number of rows
        \param ncols_b Number of columns

        If ncols_b is not given (== 0), alpha and beta have the same number of
        columns.
     **/
    ab_matrix(size_t nrows_a, size_t ncols_a,
        size_t nrows_b, size_t ncols_b = 0) : ab_object<arma::mat>(false) {

        alpha() = arma::mat(nrows_a, ncols_a);
        beta()  = arma::mat(nrows_b, ncols_b == 0 ? ncols_a : ncols_b);
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_a() const { return alpha().n_rows; }

    /** \brief Return the number of alpha-spin columns
     **/
    size_t ncols_a() const { return alpha().n_cols; }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_b() const { return beta().n_rows; }

    /** \brief Return the number of alpha-spin columns
     **/
    size_t ncols_b() const { return beta().n_cols; }

    /** \brief Add to current matrix
        \param other Matrix to add
     **/
    ab_matrix &operator+=(const ab_matrix &other) {
        alpha() += other.alpha();
        if (! is_alpha_eq_beta()) beta() += other.beta();
        return *this;
    }

    /** \brief Subtract from current matrix
        \param other Matrix to subtract
     **/
    ab_matrix &operator-=(const ab_matrix &other) {
        alpha() -= other.alpha();
        if (! is_alpha_eq_beta()) beta() -= other.beta();
        return *this;
    }

    /** \brief Scalar multiplication of the current matrix
        \param scalar Scalar factor
     **/
    ab_matrix &operator*=(double scalar) {
        alpha() *= scalar;
        if (! is_alpha_eq_beta()) beta() *= scalar;
        return *this;
    }

    /** \brief Return the transpose of the ab_matrix

        Note: this is a convenience function that does not have the full
          efficiency of directly using the armadillo operations!
     **/
    ab_matrix t() const {
        ab_matrix res(is_alpha_eq_beta());
        res.alpha() = alpha().t();
        if (! is_alpha_eq_beta()) res.beta() = beta().t();
        return res;
    }
    
    /** \brief Transpose the ab_matrix
     **/
    void inplace_trans() {
        arma::inplace_trans(alpha());
        if (! is_alpha_eq_beta()) arma::inplace_trans(beta());
    }
};

} // namespace libwfa

#endif // LIBWFA_AB_MATRIX_H
