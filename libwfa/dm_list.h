#ifndef LIBWFA_DM_LIST_H
#define LIBWFA_DM_LIST_H

#include <map>
#include "ab_matrix.h"

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class dm_type {
private:
    char m_type;

private:
    explicit dm_type(char t = 0) : m_type(t) { }

public:
    bool operator==(const dm_type &t) const { return m_type == t.m_type; }
    bool operator!=(const dm_type &t) const { return m_type != t.m_type; }

    static const dm_type sdm; //!< State density matrix type
    static const dm_type tdm; //!< Transition density matrix type
    static const dm_type adm; //!< Attachment density matrix type
    static const dm_type ddm; //!< Detachment density matrix type
    static const dm_type hdm; //!< Hole density matrix type
    static const dm_type edm; //!< Electron density matrix type
};


/** \brief Container for a set of density matrices of certain type

    The container stores the density matrix type and a list of index-density
    matrix pairs. The density matrices which are added to the list are copied
    and the copies are stored together with the unique index.

    \ingroup libwfa
 **/
class dm_list {
public:
    typedef std::map<size_t, ab_matrix *>::const_iterator iterator;

private:
    dm_type m_type; //!< Type of density matrices in the container
    std::map<size_t, ab_matrix *> m_lst; //!< List of density matrices

public:
    /** \brief Constructor
     **/
    dm_list(dm_type type) : m_type(type) { }

    /** \brief Destructor
     **/
    ~dm_list() { clear(); }

    /** \brief Add density matrix to list
        \param idx Index of density matrix
        \param dm Density matrix
     **/
    void add(size_t idx, const ab_matrix &dm);

    /** \brief Erase a density matrix
        \param idx Index of density matrix to erase
     **/
    void erase(size_t idx);

    /** \brief Empty the whole container
     **/
    void clear();

    /** \brief Rerturn the type of density store in the container
     **/
    dm_type type() const { return m_type; }

    /** \brief Number of element in list
     **/
    size_t size() const { return m_lst.size(); }

    /** \brief Return if the given index exists in the list
        \param idx Index to check
     **/
    bool index_exists(size_t idx) const {
        return m_lst.find(idx) != m_lst.end();
    }

    /** \brief Return iterator to the element with idx in list or end()
        \param idx Index
     **/
    iterator find(size_t idx) const { return m_lst.find(idx); }

    /** \brief STL-style iterator to begin of list
     **/
    iterator begin() const { return m_lst.begin(); }

    /** \brief STL-style iterator to end of list
     **/
    iterator end() const { return m_lst.end(); }

    /** \brief Index of the current density
     **/
    size_t get_index(iterator i) const { return i->first; }

    /** \brief Get current density
     **/
    const ab_matrix &get_density(iterator i) const { return *(i->second); }
};


} // namespace libwfa

#endif // LIBWFA_DM_LIST_H
