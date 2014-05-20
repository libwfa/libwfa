#ifndef LIBWFA_GRID3D_H
#define LIBWFA_GRID3D_H

namespace libwfa {


/** \brief Grid in real space (3D)

    The grid is defined by the origin and the unit vectors of one grid
    cell.

    \ingroup libwfa
 **/
class grid3d {
private:
    double m_vec[12]; //!< Origin and grid cell vectors
    unsigned int m_npts[3]; //!< # grid points in each direction

public:
    /** \brief Default constructor
        Initializes grid with
        - origin at (0,0,0)
        - cell vectors in standard direction
        - and only one point (n[i] = 1)
     **/
    grid3d();

    /** \brief Constructor
        \param n Grid points in all directions
        \param d Step size in all directions

        Initializes grid with cell vectors of length d in standard directions
        and n grid points in each direction. The grid is centered around zero.
     **/
    grid3d(unsigned int n, double d);

    /** \brief Constructor
        \param n Vector of grid points
        \param d Vector of step sizes

        Initializes grid with cell vectors of specified length (d[i]) in
        standard directions and the respective number of grid points. The
        grid is centered around (0,0,0).
     **/
    grid3d(const unsigned int (&n)[3], const double (&d)[3]);

    /** \brief Set the grid origin
        \param x0 Origin vector
     **/
    void set_origin(const double (&x0)[3]);

    /** \brief Set grid cell vector and number of points along it
        \param i Index of grid cell vector
        \param ni Number of points in the direction
        \param ei Grid cell vector
     **/
    void set_direction(unsigned int i, unsigned int ni, const double (&ei)[3]);

    /** \brief Return one coordinate of origin
        \param dim Coordinate index
     **/
    const double &origin(unsigned int dim) const;

    /** \brief Return one coordinate of one grid cell vector
        \param i Index of grid cell vector
        \param dim Coordinate index
     **/
    const double &direction(unsigned int i, unsigned int dim) const;

    /** \brief Return number of points in along one direction
        \param dir Grid direction
     **/
    unsigned int npts(unsigned int dir) const;

    /** \brief Perform consistency check
     **/
    void check() const;

private:
    static void check_idx(unsigned int dim);
};


inline const double &grid3d::origin(unsigned int dim) const {
#ifdef LIBWFA_DEBUG
    check_idx(dim);
#endif

    return m_vec[dim];
}


inline const double &grid3d::direction(unsigned int i, unsigned int dim) const {
#ifdef LIBWFA_DEBUG
    check_idx(i);
    check_idx(dim);
#endif

    return m_vec[3 + i * 3 + dim];
}

inline unsigned int grid3d::npts(unsigned int dim) const {
#ifdef LIBWFA_DEBUG
    check_idx(dim);
#endif

    return m_npts[dim];
}


} // namespace libwfa

#endif
