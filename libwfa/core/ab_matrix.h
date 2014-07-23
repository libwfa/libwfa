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
class ab_matrix : public ab_object< arma::Mat<double> > {
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta matrix
     **/
    ab_matrix(bool aeqb = false) : ab_object< arma::Mat<double> >(aeqb) { }

    /** \brief Constructor for matrix with alpha == beta
        \param nrows Number of rows
        \param ncols Number of columns

        If ncols is not given (== 0), a symmetric matrix is constructed
     **/
    ab_matrix(size_t nrows, size_t ncols = 0) :
        ab_object< arma::Mat<double> >(true) {

        alpha() = arma::Mat<double>(nrows, ncols == 0 ? nrows : ncols);
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
            size_t nrows_b, size_t ncols_b = 0) :
                ab_object< arma::Mat<double> >(false) {

        alpha() = arma::Mat<double>(nrows_a, ncols_a);
        beta()  = arma::Mat<double>(nrows_b, ncols_b == 0 ? ncols_a : ncols_b);
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

    /** \brief Add other to current matrix
     **/
    ab_matrix &operator+=(const ab_matrix &other) {
        alpha() += other.alpha();
        if (! is_alpha_eq_beta()) beta() += other.beta();
        return *this;
    }

    /** \brief Subtract other from current matrix
     **/
    ab_matrix &operator-=(const ab_matrix &other) {
        alpha() -= other.alpha();
        if (! is_alpha_eq_beta()) beta() -= other.beta();
        return *this;
    }

    /** \brief Scalar multiplication of the current matrix
     **/
    ab_matrix &operator*=(const double &scalar) {
        alpha() *= scalar;
        if (! is_alpha_eq_beta()) beta() *= scalar;
        return *this;
    }

    /** \brief Transpose the ab_matrix
     */
    ab_matrix t() const {
        ab_matrix outmat(this->is_alpha_eq_beta());

        outmat.alpha() = this->alpha().t();
        if (not this->is_alpha_eq_beta()) 
            outmat.beta() = this->beta().t();

        return outmat;
    }

    /** \brief Addition of two ab_matrix instances
     */
    ab_matrix operator+(const ab_matrix &other) const {
        bool aeqb(this->is_alpha_eq_beta() && other.is_alpha_eq_beta());
        ab_matrix outmat(aeqb);

        outmat.alpha() = this->alpha() + other.alpha();
        if (not aeqb)
            outmat.beta() = this->beta() + other.beta();

        return outmat;
    }    

    /** \brief Subtraction of two ab_matrix instances
     */
    ab_matrix operator-(const ab_matrix &other) const {
        bool aeqb(this->is_alpha_eq_beta() && other.is_alpha_eq_beta());
        ab_matrix outmat(aeqb);

        outmat.alpha() = this->alpha() - other.alpha();
        if (not aeqb)
            outmat.beta() = this->beta() - other.beta();

        return outmat;
    }

    /** \brief Matrix multiplication of two ab_matrix instances
     */
    ab_matrix operator*(const ab_matrix &other) const {
        bool aeqb(this->is_alpha_eq_beta() && other.is_alpha_eq_beta());
        ab_matrix outmat(aeqb);

        outmat.alpha() = this->alpha() * other.alpha();        
        if (not aeqb)
            outmat.beta() = this->beta() * other.beta();

        return outmat;
    }

    /** \brief Element-wise multiplication of two ab_matrix instances
     */
    ab_matrix operator%(const ab_matrix &other) const {
        bool aeqb(this->is_alpha_eq_beta() && other.is_alpha_eq_beta());
        ab_matrix outmat(aeqb);

        outmat.alpha() = this->alpha() % other.alpha();
        if (not aeqb)
            outmat.beta() = this->beta() % other.beta();

        return outmat;
    }    

};

} // namespace libwfa

#endif // LIBWFA_AB_MATRIX_H
