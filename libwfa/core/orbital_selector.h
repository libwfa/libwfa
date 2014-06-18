#ifndef LIBWFA_ORBITAL_SELECTOR_H
#define LIBWFA_ORBITAL_SELECTOR_H

#include "selector.h"

namespace libwfa {

/** \brief Selects orbital indexes in the range [0,n] as occupied and virtual
        orbitals.

    \ingroup libwfa
 **/
class orbital_selector {
private:
    selector m_occ, m_vir; //!< Selectors for occupied and virtual indexes

public:
    /** \brief Constructor
        \param ntotal Number of elements in range

        Range starts at 0 and selection is empty.
     **/
    orbital_selector(size_t ntotal = 0) : m_occ(ntotal), m_vir(ntotal) { }

    /** \brief Add one index to selection
     **/
    void select(bool as_occ, size_t i);

    /** \brief Add range of indexes to selection
        \param as_occ Are the indexes occupied indexes
        \param begin Start \f$ i_0 \f$ of first range of selected elements
        \param end End \f$ i_1 \f$ of first range of selected elements
        \param inc Increment \f$ \Delta \f$ in range
        \param reverse If true, go from end to begin and select elements
     **/
    void select(bool as_occ, size_t begin, size_t end, size_t inc = 1,
            bool reverse = false);

    /** \brief Deselect one index in selection
     **/
    void deselect(size_t i);

    /** \brief Return number of elements in range
     **/
    size_t n_indexes() const { return m_occ.n_indexes(); }

    /** \brief Return number of selected elements
     **/
    size_t n_selected(bool occ) const {
        return (occ ? m_occ.n_selected() : m_vir.n_selected());
    }

    /** \brief Return number of selected elements
     **/
    size_t n_selected() const {
        return m_occ.n_selected() + m_vir.n_selected();
    }

    /** \brief Are all elements selected?
     **/
    bool all_selected() const {
        return m_occ.n_selected() + m_vir.n_selected() == m_occ.n_indexes();
    }

    /** \brief Are no elements selected?
     **/
    bool none_selected() const {
        return m_occ.none_selected() && m_vir.none_selected();
    }

    /** \brief Is i selected
     **/
    bool is_selected(bool as_occ, size_t i) const {
        return (as_occ ? m_occ.is_selected(i) : m_vir.is_selected(i));
    }

    /** \brief Is i selected
     **/
    bool is_selected(size_t i) const {
        return m_occ.is_selected(i) || m_vir.is_selected(i);
    }

    /** \brief Construct a vector of selected elements
     **/
    std::vector<size_t> get_selected() const;

    /** \brief Construct a vector of selected elements (return Armadillo vector)
     **/
    arma::Col<arma::uword> get_selected_arma() const;
};

} // namespace libwfa

#endif // LIBWFA_ORBITAL_SELECTOR_H
