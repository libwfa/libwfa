#ifndef LIBWFA_CUBE_WRITER_H
#define LIBWFA_CUBE_WRITER_H

#include <fstream>
#include "atom_list.h"
#include "grid3d.h"

namespace libwfa {


/** \brief Write data as cube file

    \ingroup libwfa
 **/
class cube_writer {
private:
    std::ofstream m_out; //!< Output file
    size_t m_npoints; //!< Number of volumetric data points
    size_t m_nrec; //!< Length of one data record (usually nz of grid)
    size_t m_cur; //!< Current point in volumetric data

public:
    /** \brief Constructor
        \param filename File name
        \param line1 First comment line
        \param line2 Second comment line
        \param grid Grid data
        \param atoms List of atoms
     **/
    cube_writer(const std::string &filename,
        const std::string &line1, const std::string &line2,
        const grid3d &grid, const atom_list &atoms);

    /** \brief Return the total number of points to be written
     **/
    size_t npoints() const { return m_npoints; }

    /** \brief Return if writing of cube file is complete
     **/
    bool complete() const { return m_cur == m_npoints; }

    /** \brief Write batch of data to disk
        \param nbatch Size of data
        \param data Data array
        \return Number of points written
     */
    size_t write(size_t nbatch, const double *data);

private:
    void write_header(const std::string &line1, const std::string &line2,
        const grid3d &grid, const atom_list &atoms);

};


}

#endif // LIBWFA_CUBE_WRITER_H
