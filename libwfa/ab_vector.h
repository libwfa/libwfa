#ifndef LIBWFA_AB_VECTOR_H
#define LIBWFA_AB_VECTOR_H

#include <armadillo>

namespace libwfa {

/** \brief Container for spin-vectors

    \ingroup libwfa
 **/
class ab_vector {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    arma::Col<double> m_vec_a; //!< Alpha-spin vector
    arma::Col<double> m_vec_b; //!< Beta-spin vector
    
public:
    /** \brief Constructor
        \param nrows_a Number of alpha-spin rows
        \param nrows_b Number of beta-spin rows
        \param aeqb True, if alpha and beta spin vector are
     **/
    ab_vector(size_t nrows_a = 0, size_t nrows_b = 0, bool aeqb = false) :
        m_vec_a(nrows_a), m_vec_b(aeqb ? 0 : nrows_b), m_aeqb(aeqb) {
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_a() const {
        return m_vec_a.n_rows;
    }

    /** \brief Return the number of beta-spin rows
     **/
    size_t nrows_b() const {
        return (m_aeqb ? m_vec_a.n_rows : m_vec_b.n_rows);
    }

    /** \brief Change if alpha == beta
     **/
    void set_aeqb(bool aeqb) {
        if (aeqb == m_aeqb) return;

        if (m_aeqb) { m_vec_b = m_vec_a; }
        else { m_vec_b.resize(0, 0); }
        m_aeqb = aeqb;
    }

    /** \brief Are alpha- and beta-spin matrices identical
     **/
    bool is_alpha_eq_beta() const {
        return m_aeqb;
    }

    /** \brief Return alpha-spin matrix
     **/
    arma::Col<double> &alpha() {
        return m_vec_a;
    }
    
    /** \brief Return alpha-spin matrix (const version)
     **/
    const arma::Col<double> &alpha() const {
        return m_vec_a;
    }
    
    /** \brief Return beta-spin matrix
     **/
    arma::Col<double> &beta() {
        return (m_aeqb ? m_vec_a : m_vec_b);
    }
    
    /** \brief Return beta-spin matrix (const version
     **/
    const arma::Col<double> &beta() const {
        return (m_aeqb ? m_vec_a : m_vec_b);
    }
    
};

} // namespace adcman

#endif // LIBWFA_AB_VECTOR_H
