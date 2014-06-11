#ifndef LIBWFA_SELECTOR_H
#define LIBWFA_SELECTOR_H

#include <armadillo>
#include <cstdlib>
#include <map>

namespace libwfa {

/** \brief Selects indexes in the range [0,n]

    \ingroup libwfa
 **/
class selector {
private:
    std::map<size_t, size_t> m_indexes; //!< Selected elements (
    size_t m_ntotal; //!< Total number of indexes
    size_t m_cur; //!< Current index

public:
    /** \brief Constructor
        \param ntotal Number of elements in range

        Range starts at 0 and selection is empty.
     **/
    selector(size_t ntotal = 0) : m_ntotal(ntotal), m_cur(0) { }

    /** \brief Add one index to selection
     **/
    void select(size_t i);

    /** \brief Add range of indexes to selection
        \param begin Start \f$ i_0 \f$ of first range of selected elements
        \param end End \f$ i_1 \f$ of first range of selected elements
        \param inc Increment \f$ \Delta \f$ in range
        \param reverse If true, go from end to begin and select elements
     **/
    void select(size_t begin, size_t end, size_t inc = 1, bool reverse = false);

    /** \brief Select all elements (in the given order)
     **/
    void select_all();

    /** \brief Deselect one index in selection
     **/
    void deselect(size_t i);

    /** \brief Deselect all indexes
     **/
    void deselect_all();

    /** \brief Return number of elements in range
     **/
    size_t n_indexes() const { return m_ntotal; }

    /** \brief Return number of selected elements
     **/
    size_t n_selected() const { return m_indexes.size(); }

    /** \brief Are all elements selected?
     **/
    bool all_selected() const { return m_ntotal == m_indexes.size(); }

    /** \brief Are no elements selected?
     **/
    bool none_selected() const { return m_indexes.size() == 0; }

    /** \brief Is i selected
     **/
    bool is_selected(size_t i) const {
        check(i);
        return m_indexes.count(i) != 0;
    }

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
