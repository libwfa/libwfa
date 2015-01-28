#ifndef LIBWFA_GRID3D_H
#define LIBWFA_GRID3D_H

#include <armadillo>
#include <cstddef>

namespace libwfa {


/** \brief Grid in real space (3D)

    The grid is defined by the origin and the unit vectors of one grid
    cell.

    \ingroup libwfa
 **/
class grid3d {
private:
    arma::mat m_vec; //!< Origin and direction (column) vectors (4x3)
    arma::Col<unsigned int> m_npts; //!< Column vector with points
    size_t m_ntotal; //!< Total number of points

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
    grid3d(const arma::Col<unsigned int> &n, const arma::vec &d);

    /** \brief Set the grid origin
        \param x0 Origin vector
     **/
    void set_origin(const arma::vec &x0);

    /** \brief Set grid cell vector and number of points along it
        \param i Index of grid cell vector
        \param ni Number of points in the direction
        \param ei Grid cell vector
     **/
    void set_direction(unsigned int i,
            unsigned int ni, const arma::vec &ei);

    /** \brief Returns the total number of grid points
     **/
    size_t size() const { return m_ntotal; }

    /** \brief Return origin column vector
     **/
    arma::subview<double> origin() const {
        return m_vec.col(0);
    }

    /** \brief Return one grid cell vector
        \param i Index of grid cell vector
     **/
    arma::subview<double> direction(unsigned int i) const;

    /** \brief Return number of points in along one direction
     **/
    const arma::Col<unsigned int> &npts() const {
        return m_npts;
    }

    /** \brief Perform consistency check
     **/
    void check() const;

    /** \brief Build 3D coordinates for a set of grid points
        \param[in] i0 Index of first grid point to compute
        \param[out] pts Data array to store points (expected size: Mat(3, sz))
        \return Number of points built

        Computes sz grid points starting with the i0-th grid point. The
        algorithm assumes the usual linearization of 3d indexes using
        the last direction as the running index, i.e.
        \f$ (i,j,k) \rightarrow (i\cdot n_1 + j) \cdot n_2 + k \f$

        Does not resize the array!
     **/
    size_t build_pts(size_t i0, arma::mat &pts) const;

private:
    static void check_idx(unsigned int dim);
};


inline arma::subview<double> grid3d::direction(unsigned int i) const {
#ifdef LIBWFA_DEBUG
    check_idx(i);
#endif

    return m_vec.col(1 + i);
}


} // namespace libwfa

#endif
