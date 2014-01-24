#ifndef LIBWFA_MOLDEN_FILE_BASE_H
#define LIBWFA_MOLDEN_FILE_BASE_H

#include "ab_matrix.h"
#include "ab_vector.h"

namespace libwfa {


/** \brief Grid in real space (3D)

    \ingroup libwfa
 **/
struct grid3d {
    double ranges[6]; //!< [min, max] ranges in each direction
    unsigned int npts[3]; //!< Number of grid points in each direction

    /** \brief Constructor
     **/
    grid3d();

    //! \name Setters (with error checks)
    //@{
    void set_xrange(double min, double max, unsigned int npts);
    void set_yrange(double min, double max, unsigned int npts);
    void set_zrange(double min, double max, unsigned int npts);
    //@}

    /** \brief Consistency check
     **/
    void check() const;
};


/** \brief Base class to export data as cube files.

    Base class to export certain data (see data_type for possible data) as
    grid data into cube files.

    \ingroup libwfa
 **/
struct export_cube_base {

    grid3d grid; //!< Grid to use

    /** \brief Constructor
        \param g Grid of the cube
     **/
    export_cube_base(const grid3d &g) : grid(g) { }

    /** \brief Virtual destructor
     **/
    virtual ~export_cube_base() { }

    /** \brief Evaluate an matrix in AO basis on the grid and export as cube
        \param name Name associated with the matrix
        \param mat Matrix data in AO basis
     **/
    virtual void perform(const std::string &name,
        const arma::Mat<double> &mat) = 0;

    /** \brief Evaluate a set of vectors in AO basis on the grid and export as cube
        \param prefix Name associated with the vectors (used as prefix)
        \param idx Vector of indexes
        \param vecs Set of vectors in AO basis
     **/
    virtual void perform(const std::string &prefix,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs) = 0;

};

} // namespace adcman

#endif
