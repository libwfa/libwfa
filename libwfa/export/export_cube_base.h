#ifndef LIBWFA_EXPORT_CUBE_BASE_H
#define LIBWFA_EXPORT_CUBE_BASE_H

#include <map>
#include "export_cube_i.h"
#include "grid3d.h"

namespace libwfa {


/** \brief Base class to export data as Gaussian cube files.

    Base class to export density matrices and orbitals expressed in terms of
    basis functions as grid data into cube files. Any derived class has to
    implement the method \c evaluate_on_grid(...) (see below).

    The class implements the export by first collecting the data sets and
    exporting them only when the number of data sets exceeds \c m_nmax or the
    function \c do_export() is called.

    The grid data of the density matrices and orbitals are exported into
    separate files, one per density matrix or orbital. The naming scheme of
    the files is as follows: for density matrices the file name is constructed
    as
    \code prefix + name + ".cube" \endcode
    while for orbitals the index is appended as well
    \code prefix + name + index + ".cube" \endcode

    \ingroup libwfa
 **/
class export_cube_base : public export_cube_i {
public:
    static const char k_clazz[]; //!< Class name

protected:
    //! \brief Variables to set by child class
    //@{
    std::string m_comment; //!< Comment
    size_t m_batchsz; //!< Batch size
    size_t m_nmax; //!< Max number of data sets to be stored
    //@}

private:
    struct dm_data {
        std::string desc; //!< Description
        arma::mat data; //!< Density matrix data

        dm_data(const std::string &d, const arma::mat &dm) :
            desc(d), data(dm) { }
    };
    struct orb_data {
        std::string desc; //!< Description
        std::vector<size_t> idx; //!< Indexes
        arma::mat data; //!< Orbital data

        orb_data(const std::string &d, const std::vector<size_t> &i,
            const arma::mat &orb) : desc(d), idx(i), data(orb) { }
    };
    typedef std::map<std::string, dm_data *> dm_list;
    typedef std::map<std::string, orb_data *> orb_list;

private:
    grid3d m_grid; //!< The grid
    const arma::uvec &m_atnum; //!< Atom numbers (dim: N)
    cosnt arma::mat &m_coords; //!< Atomic coordinates (dim: 3 x N)
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
    export_cube_base(const grid3d &grid, const arma::uvec &atnum,
        const arma::mat &coord, const std::string prefix = "");

    /** \brief Destructor
     **/
    virtual ~export_cube_base() {
        clear_data();
    }

    /** \brief Return the grid data
     **/
    const grid3d &grid() const { return m_grid; }

    /** \brief Return the vector of atomic numbers
     **/
    const arma::uvec &atomic_numbers() const { return m_atnum; }

    /** \brief Return the 3xN matrix of atom positions
     **/
    const arma::mat &atomic_coordinates() const { return m_coords; }

    /** \copydoc export_cube_i::perform
     **/
    virtual void perform(const std::string &name, const std::string &desc,
        const arma::mat &mat);

    /** \copydoc export_cube_i::perform(const std::string&, const std::string&, const std::vector<size_t>&, const arma::mat&)
     **/
    virtual void perform(const std::string &name, const std::string &desc,
        const std::vector<size_t> &idx, const arma::mat &vecs);

    /** \brief Export the stored data as cube files
     **/
    void do_export();

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
    virtual void evaluate_on_grid(const arma::mat &pts, arma::mat &b2g) = 0;

private:
    void clear_data();
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_CUBE_BASE_H
