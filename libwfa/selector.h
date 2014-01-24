#ifndef LIBWFA_SELECTOR_H
#define LIBWFA_SELECTOR_H

#include <armadillo>
#include <cstdlib>
#include <vector>

namespace libwfa {

/** \brief Selects a number of indexes in the range [0,n]

    \ingroup libwfa
 **/
class selector {
private:
    std::vector<bool> m_indexes; //!< Selected/unselected flag for each element
    size_t m_nselected; //!< Number of selected elements

public:
    /** \brief Constructor
        \param ntotal Number of elements in range

        Range starts at 0 and selection is empty.
     **/
    selector(size_t ntotal = 0) : m_indexes(ntotal, false), m_nselected(0) { }

    /** \brief Add one index to selection
     **/
    void select(size_t i);

    /** \brief Add range of indexes to selection
        \param begin Start \f$ i_0 \f$ of first range of selected elements
        \param end End \f$ i_1 \f$ of first range of selected elements
        \param inc Increment \f$ \Delta \f$ in range
     **/
    void select(size_t begin, size_t end, size_t inc = 1);

    /** \brief Select all elements
     **/
    void select_all();

    /** \brief Deselect one index in selection
     **/
    void deselect(size_t i);

    /** \brief Deselect one index in selection
     **/
    void deselect_all();

    /** \brief Return number of elements in range
     **/
    size_t n_indexes() const { return m_indexes.size(); }

    /** \brief Return number of selected elements
     **/
    size_t n_selected() const { return m_nselected; }

    /** \brief Are all elements selected?
     **/
    bool all_selected() const { return m_nselected == m_indexes.size(); }

    /** \brief Are no elements selected?
     **/
    bool none_selected() const { return m_nselected == 0; }

    /** \brief Is i selected
     **/
    bool is_selected(size_t i) const { check(i); return m_indexes[i]; }

    /** \brief Construct a vector of selected elements
     **/
    std::vector<size_t> get_selected() const;

    /** \brief Construct a vector of selected elements (return Armadillo vector)
     **/
    arma::Col<arma::uword> get_selected_arma() const;

private:
    void check(size_t i) const;

};

} // namespace libwfa

#endif // LIBWFA_SELECTOR_H
