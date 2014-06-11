#ifndef LIBWFA_EXPORT_CUBE_H
#define LIBWFA_EXPORT_CUBE_H

#include <armadillo>
#include <map>
#include <string>
#include "grid3d.h"

namespace libwfa {


/** \brief Base class to export data as Gaussian cube files.

    Base class to export density matrices and orbitals expressed in terms of
    basis functions as grid data into cube files. Any derived classes have to
    implement the method \c evaluate_on_grid(...) (see below).

    The grid data of the density matrices and orbitals is exported into
    separate files, one per density matrix or orbital. The naming scheme of
    the files is as follows: for density matrices the file name is constructed
    as
    \code prefix + name + ".cube" \endcode
    while for orbitals the index is appended as well
    \code prefix + name + index + ".cube" \endcode

    \ingroup libwfa
 **/
class export_cube_base {
public:
    static const char k_clazz[]; //!< Class name

protected:
    //! \brief Variables to set by child class
    //@{
    std::string m_comment; //!< Comment
    size_t m_batchsz; //!< Batch size
    //@}

private:
    struct dm_data {
        std::string desc; //!< Description
        arma::Mat<double> data; //!< Density matrix data

        dm_data(const std::string &d, const arma::Mat<double> &dm) :
            desc(d), data(dm) { }
    };
    struct orb_data {
        std::string desc; //!< Description
        std::vector<size_t> idx; //!< Indexes
        arma::Mat<double> data; //!< Orbital data

        orb_data(const std::string &d, const std::vector<size_t> &i,
            const arma::Mat<double> &orb) : desc(d), idx(i), data(orb) { }
    };
    typedef std::map<std::string, dm_data *> dm_list;
    typedef std::map<std::string, orb_data *> orb_list;

private:
    const grid3d &m_grid; //!< The grid
    const arma::Col<unsigned int> &m_atnum; //!< Atom numbers (dim: N)
    const arma::Mat<double> &m_coords; //!< Atomic coordinates (dim: 3 x N)
    dm_list m_dms; //!< Density matrices to export as cube
    orb_list m_orbs; //!< Orbitals to export as cube
    std::string m_prefix; //!< Filename prefix

public:
    /** \brief Constructor
        \param grid Grid to generate the volumetric data
        \param atnum List of atomic numbers
        \param coords List of atomic coords (#atoms columns x 3 rows)
        \param prefix Prefix to use for the filenames (e.g. directory)
     **/
    export_cube_base(const grid3d &grid, const arma::Col<unsigned int> &atnum,
        const arma::Mat<double> &coord, const std::string prefix = "");

    /** \brief Destructor
     **/
    virtual ~export_cube_base() {
        clear_data();
    }

    /** \brief Add a density matrix in terms of basis functions for evaluation
            on a grid
        \param name Name associated with the matrix (used as filename)
        \param desc Description of the matrix (to put as comment, one line)
        \param mat Density matrix in terms of basis functions
     **/
    void add(const std::string &name, const std::string &desc,
        const arma::Mat<double> &mat);

    /** \brief Add a set of vectors in terms of basis functions for evaluation
            on a grid
        \param name Name associated with the vectors (used as filename)
        \param desc Description of the vectors (to put as comment, one line)
        \param vecs Set of vectors in basis functions (column vectors)
     **/
    void add(const std::string &name, const std::string &desc,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs);

    /** \brief Export the stored data as cube files
     **/
    void perform();

protected:
    /** \brief Evaluate basis functions on a number of grid points
        \param[in] pts Grid points (NPT columns of length 3)
        \param[out] b2g Values of atomic basis functions on grid

        Implementations of the function are expected to return a matrix
        with NAO columns each of length NPT, where NAO is the number of atomic
        basis functions and NPT is the number of grid points. Each column
        should contain the respective basis function evaluated at the grid
        points provided.
     **/
    virtual void evaluate_on_grid(const arma::Mat<double> &pts,
            arma::Mat<double> &b2g) = 0;

private:
    void clear_data();
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_CUBE_H
