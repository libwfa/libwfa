#ifndef LIBWFA_AB_MATRIX_H
#define LIBWFA_AB_MATRIX_H

#include <armadillo>

namespace libwfa {

/** \brief Container for spin-matrices

    \ingroup libwfa
 **/
class ab_matrix {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    arma::Mat<double> m_mat_a; //!< Alpha-spin matrix
    arma::Mat<double> m_mat_b; //!< Beta-spin matrix
    
public:
    /** \brief Constructor
        \param nrows_a Number of rows
        \param ncols_a Number of columns
        \param nrows_b Number of rows
        \param ncols_b Number of columns
        \param aeqb True, if alpha and beta spin matrices are
     **/
    ab_matrix(
        size_t nrows_a = 0, size_t ncols_a = 0,
        size_t nrows_b = 0, size_t ncols_b = 0,
        bool aeqb = false) :
        m_aeqb(aeqb && nrows_a == nrows_b && ncols_a == ncols_b),
        m_mat_a(nrows_a, ncols_a),
        m_mat_b(m_aeqb ? 0 : nrows_b, m_aeqb ? 0 : ncols_b)
         {
    }

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
        return (m_aeqb ? m_mat_a.n_cols : m_mat_b.n_rows);
    }

    /** \brief Change if alpha == beta
     **/
    void set_aeqb(bool aeqb) {
        if (aeqb == m_aeqb) return;

        if (m_aeqb) { m_mat_b = m_mat_a; }
        else { m_mat_b.resize(0, 0); }
        m_aeqb = aeqb;
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
