#ifndef LIBWFA_EXPORT_CUBE_QCHEM_H
#define LIBWFA_EXPORT_CUBE_QCHEM_H

#include "../export_cube_i.h"

namespace libwfa {


/** \brief Base class to export data as cube files.

    Base class to export certain data (see data_type for possible data) as
    grid data into cube files.

    \ingroup libwfa
 **/
class export_cube_qchem : public export_cube_i {
private:
    const grid3d &m_grid; //!< Grid to use
    const atom_list &m_atoms; //!< Atom lists
    std::string m_path; //!< Path where to put files

public:
    /** \brief Constructor
        \param grid Grid to generate the volumetric data
        \param atoms List of atoms
        \param path Path where to put the cube files
     **/
    export_cube_qchem(const grid3d &grid, const atom_list &atoms,
        const std::string path = "") :
        m_grid(grid), m_atoms(atoms), m_path(path)
    { }

    /** \brief Virtual destructor
     **/
    virtual ~export_cube_qchem() { }

    /** \brief Evaluate an matrix in AO basis on the grid and export as cube
        \param name Name associated with the matrix
        \param mat Matrix data in AO basis
     **/
    virtual void perform(const std::string &name,
        const arma::Mat<double> &mat);

    /** \brief Evaluate a set of vectors in AO basis on the grid and export
            as cube
        \param prefix Name associated with the vectors (used as prefix)
        \param idx Vector of indexes
        \param vecs Set of vectors in AO basis
     **/
    virtual void perform(const std::string &prefix,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs);
private:
    void build_grid_points(double *gpts, size_t start, size_t sz)
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_CUBE_QCHEM_H
