#ifndef LIBWFA_AB_SELECTOR_H
#define LIBWFA_AB_SELECTOR_H

#include "selector.h"

namespace libwfa {

/** \brief Container for spin-selectors

    \ingroup libwfa
 **/
class ab_selector {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    selector m_sel_a; //!< Alpha-spin selector
    selector m_sel_b; //!< Beta-spin selector
    
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta spin vector
     **/
    ab_selector(bool aeqb = false) : m_aeqb(aeqb) { }


    /** \brief Constructor for alpha == beta
        \param nidx Number of indexes in range
     **/
    ab_selector(size_t nidx) :
        m_aeqb(true), m_sel_a(nidx), m_sel_b(0) { }

    /** \brief Constructor for alpha != beta
        \param nidx_a Number of alpha-spin indexes
        \param nidx_b Number of beta-spin indexes
     **/
    ab_selector(size_t nidx_a, size_t nidx_b) :
        m_aeqb(false), m_sel_a(nidx_a), m_sel_b(nidx_b) {
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nidx_a() const {
        return m_sel_a.n_indexes();
    }

    /** \brief Return the number of beta-spin rows
     **/
    size_t nidx_b() const {
        return (m_aeqb ? m_sel_a.n_indexes() : m_sel_b.n_indexes());
    }

    /** \brief Set alpha == beta
     **/
    void set_alpha_eq_beta() {
        if (m_aeqb) return;
        m_sel_b = selector(0);
        m_aeqb = true;
    }

    /** \brief Set alpha != beta
     **/
    void set_alpha_neq_beta() {
        if (! m_aeqb) return;
        m_sel_b = m_sel_a;
        m_aeqb = false;
    }

    /** \brief Are alpha- and beta-spin matrices identical
     **/
    bool is_alpha_eq_beta() const {
        return m_aeqb;
    }

    /** \brief Return alpha-spin matrix
     **/
    selector &alpha() {
        return m_sel_a;
    }
    
    /** \brief Return alpha-spin matrix (const version)
     **/
    const selector &alpha() const {
        return m_sel_a;
    }
    
    /** \brief Return beta-spin matrix
     **/
    selector &beta() {
        return (m_aeqb ? m_sel_a : m_sel_b);
    }
    
    /** \brief Return beta-spin matrix (const version
     **/
    const selector &beta() const {
        return (m_aeqb ? m_sel_a : m_sel_b);
    }
    
};

} // namespace adcman

#endif // LIBWFA_AB_SELECTOR_H
