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
    /** \brief Default constructor
        \param aeqb If true, alpha == beta spin vector
     **/
    ab_vector(bool aeqb = false) : m_aeqb(aeqb) { }


    /** \brief Constructor for alpha == beta
        \param nrows Number of rows
     **/
    ab_vector(size_t nrows) :
        m_aeqb(true), m_vec_a(nrows) { }

    /** \brief Constructor for alpha != beta
        \param nrows_a Number of alpha-spin rows
        \param nrows_b Number of beta-spin rows
     **/
    ab_vector(size_t nrows_a, size_t nrows_b) :
        m_aeqb(false), m_vec_a(nrows_a), m_vec_b(nrows_b) {
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

    /** \brief Set alpha == beta
     **/
    void set_alpha_eq_beta() {
        if (m_aeqb) return;
        m_vec_b.resize(0, 0);
        m_aeqb = true;
    }

    /** \brief Set alpha != beta
     **/
    void set_alpha_neq_beta() {
        if (! m_aeqb) return;
        m_vec_b = m_vec_a;
        m_aeqb = false;
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
