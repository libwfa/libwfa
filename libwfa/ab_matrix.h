#ifndef LIBWFA_AB_MATRIX_H
#define LIBWFA_AB_MATRIX_H

#include <armadillo>

namespace libwfa {

/** \brief Container for matrices with spin property

    The container stores two matrices with alpha- and beta-spin and a flag to
    indicate, if both matrices have to be identical (i.e. spin-restricted
    calculations). The case of alpha == beta will be enforced by the container
    by using the same matrix for alpha-spin and beta-spin.

    \ingroup libwfa
 **/
class ab_matrix {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    arma::Mat<double> m_mat_a; //!< Alpha-spin matrix
    arma::Mat<double> m_mat_b; //!< Beta-spin matrix
    
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta matrix
     **/
    ab_matrix(bool aeqb = false) : m_aeqb(aeqb) { }

    /** \brief Constructor for matrix with alpha == beta
        \param nrows Number of rows
        \param ncols Number of columns

        If ncols is not given (== 0), a symmetric matrix is constructed
     **/
    ab_matrix(size_t nrows, size_t ncols = 0) :
        m_aeqb(true), m_mat_a(nrows, (ncols == 0 ? nrows : ncols)) { }

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
        m_aeqb(false), m_mat_a(nrows_a, ncols_a),
        m_mat_b(nrows_b, (ncols_b == 0 ? ncols_a : ncols_b)) { }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_a() const {
        return m_mat_a.n_rows;
    }

    /** \brief Return the number of alpha-spin columns
     **/
    size_t ncols_a() const {
        return m_mat_a.n_cols;
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_b() const {
        return (m_aeqb ? m_mat_a.n_rows : m_mat_b.n_rows);
    }

    /** \brief Return the number of alpha-spin columns
     **/
    size_t ncols_b() const {
        return (m_aeqb ? m_mat_a.n_cols : m_mat_b.n_cols);
    }

    /** \brief Set alpha == beta

        Modifies the container to enforce both matrices to be identical. The
        contents of the beta matrix is deleted and only the alpha matrix is
        kept.
     **/
    void set_alpha_eq_beta() {
        if (m_aeqb) return;
        m_mat_b.resize(0, 0);
        m_aeqb = true;
    }

    /** \brief Set alpha != beta

        Modifies the container to allow the matrices with alpha and beta-spin
        to change independently. However, both matrices will be identical
        copies of one another first.
     **/
    void set_alpha_neq_beta() {
        if (! m_aeqb) return;
        m_mat_b = m_mat_a;
        m_aeqb = false;
    }

    /** \brief Are alpha- and beta-spin matrices identical
     **/
    bool is_alpha_eq_beta() const {
        return m_aeqb;
    }

    /** \brief Return alpha-spin matrix
     **/
    arma::Mat<double> &alpha() {
        return m_mat_a;
    }
    
    /** \brief Return alpha-spin matrix (const version)
     **/
    const arma::Mat<double> &alpha() const {
        return m_mat_a;
    }
    
    /** \brief Return beta-spin matrix
     **/
    arma::Mat<double> &beta() {
        return (m_aeqb ? m_mat_a : m_mat_b);
    }
    
    /** \brief Return beta-spin matrix (const version
     **/
    const arma::Mat<double> &beta() const {
        return (m_aeqb ? m_mat_a : m_mat_b);
    }
    
};

} // namespace adcman

#endif // LIBWFA_AB_MATRIX_H
