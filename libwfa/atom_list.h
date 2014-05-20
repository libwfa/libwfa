#ifndef LIBWFA_ATOM_LIST_H
#define LIBWFA_ATOM_LIST_H

#include <iostream>

namespace libwfa {


/** \brief Simple wrapper around atomic coordinates

    Do not copy the objects

    \ingroup libwfa
 **/
class atom_list {
private:
    size_t m_natoms; //!< Number of atoms
    int *m_atnum; //!< Atomic numbers
    double *m_coords; //!< Atomic coordinates
    bool m_alloc; //!< Has the data has been allocated

public:
    /** \brief Constructor
        \param natoms Number of atoms
     **/
    atom_list(size_t natoms);

    /** \brief Constructor copying the data
        \param natoms Number of atoms
        \param atnum Atomic numbers
        \param coords Atomic coordinates
     **/
    atom_list(size_t natoms, const int *atnum, const double *coords);

    /** \brief Constructor
        \param natoms Number of atoms
        \param atnum Atomic numbers
        \param coords Atomic coordinates
        \param do_copy Copy the data
     **/
    atom_list(size_t natoms, int *atnum, double *coords, bool do_copy);

    /** \brief Destructor
     **/
    ~atom_list();

    /** \brief Set the atom
        \param idx Index of atom
        \param atnum Atomic number
        \param x Atomic coordinate
     **/
    void set_atom(size_t idx, int atnum, const double (&x)[3]);

    /** \brief Return the number of atoms
     **/
    size_t natoms() const { return m_natoms; }

    /** \brief Retrieve atomic number
        \param idx Index of atom
     **/
    int atomic_number(size_t idx) const;

    /** \brief Coordinate of atom
        \param idx Index of atom
        \param dim Coordinate
     **/
    const double &position(size_t idx, unsigned int dim) const;

    /** \brief Position vector of atom
        \param idx Index of atom
     **/
    const double *position(size_t idx) const;

private:
    static void check_idx(size_t idx, size_t limit);
};


/** \brief Print atom list
    \param out Output stream
    \param lst Atom list
    \return Output stream

    Formatted print in 4 columns:
    atomic_number x-coord y-coord z-coord
 **/
std::ostream &operator<<(std::ostream &out, const atom_list &lst);


inline int atom_list::atomic_number(size_t idx) const {
#ifdef LIBWFA_DEBUG
    check_idx(idx, m_natoms);
#endif
    return m_atnum[idx];
}


inline const double &atom_list::position(size_t idx, unsigned int dim) const {
#ifdef LIBWFA_DEBUG
    check_idx(idx, m_natoms);
    check_idx(dim, 3);
#endif
    return m_coords[idx * 3 + dim];
}


inline const double *atom_list::position(size_t idx) const {
#ifdef LIBWFA_DEBUG
    check_idx(idx, m_natoms);
#endif
    return m_coords + idx * 3;
}


} // namespace libwfa

#endif // LIBWFA_ATOM_LIST_H
