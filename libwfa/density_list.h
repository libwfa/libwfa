#ifndef LIBWFA_DENSITY_LIST_H
#define LIBWFA_DENSITY_LIST_H

#include <list>
#include "ab_matrix.h"

namespace libwfa {


/** \brief Enum class of density types

    \ingroup libwfa
 **/
class density_type {
private:
    char m_type;

private:
    explicit density_type(char t = 0) : m_type(t) { }

public:
    bool operator==(const density_type &t) const { return m_type == t.m_type; }
    bool operator!=(const density_type &t) const { return m_type != t.m_type; }

    static const density_type state_density; //!< State density type
    static const density_type trans_density; //!< Transition density type
    static const density_type att_density; //!< Attachment density type
    static const density_type det_density; //!< Detachment density type
};


/** \brief Container for a list of density matrices

    The container does not create copies of the density matrices, but stores
    the references. Thus, the user should make sure that the matrix objects
    will not be destroyed before the container is.

    \ingroup libwfa
 **/
class density_list {
private:
    struct dm_container {
        size_t idx; //!< Index of the density
        const ab_matrix &dm; //!< Density matrix

        dm_container(size_t idx_, const ab_matrix &dm_) :
            idx(idx_), dm(dm_) { }
    };

public:
    typedef std::list<dm_container *>::const_iterator iterator;

private:
    density_type m_type; //!< Type of densities in the container
    std::list<dm_container *> m_lst; //!< List of densities

public:
    /** \brief Constructor
     **/
    density_list(density_type type) : m_type(type) { }

    /** \brief Destructor
     **/
    ~density_list();

    /** \brief Add density matrix to list
        \param idx Index of density matrix
        \param dm Density matrix
     **/
    void add(size_t idx, const ab_matrix &dm);

    /** \brief Rerturn the type of density store in the container
     **/
    density_type type() const { return m_type; }

    /** \brief Number of element in list
     **/
    size_t size() const { return m_lst.size(); }

    /** \brief STL-style iterator to begin of list
     **/
    iterator begin() const { return m_lst.begin(); }

    /** \brief STL-style iterator to end of list
     **/
    iterator end() const { return m_lst.end(); }

    /** \brief Index of the current density
     **/
    size_t get_index(iterator i) const { return (*i)->idx; }

    /** \brief Get current density
     **/
    const ab_matrix &get_density(iterator i) const { return (*i)->dm; }

};

} // namespace libwfa


#endif // LIBWFA_DENSITY_LIST_H
