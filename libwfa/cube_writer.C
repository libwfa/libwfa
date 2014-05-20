#include <iomanip>
#include "cube_writer.h"
#include "libwfa_exception.h"

namespace libwfa {


cube_writer::cube_writer(const std::string &filename,
    const std::string &line1, const std::string &line2,
    const grid3d &grid, const atom_list &atoms) :
    m_out(filename), m_npoints(1), m_nrec(grid.npts(2)), m_cur(0) {

    grid.check();

    for (unsigned int i = 0; i < 3; i++) m_npoints *= grid.npts(i);

    write_header(line1, line2, grid, atoms);
}


size_t cube_writer::write(size_t nbatch, const double *data) {

    // Set the correct output format for data
    m_out << std::setw(13) << std::setprecision(5) << std::scientific;

    const double *ptr = data;
    size_t i = 0, j = m_cur % m_nrec, k = j % 6;
    while (i < nbatch && m_cur < m_npoints) {
        m_out << " " << ptr++;
        i++; m_cur++;
        j++; k++;

        if (j == m_nrec) {
            m_out << std::endl;
            j = 0; k = 0;
        }
        else if (k == 6) {
            m_out << std::endl;
            k = 0;
        }
    }

    return i;
}


void cube_writer::write_header(const std::string &l1,
    const std::string &l2, const grid3d &grid, const atom_list &atoms) {

    m_out << l1 << std::endl;
    m_out << l2 << std::endl;
    m_out << std::setw(5) << atoms.natoms();
    m_out << std::setw(12) << std::setprecision(6) << std::fixed;
    for (size_t i = 0; i < 3; i++) m_out << grid.origin(i);
    m_out << std::endl;
    for (size_t i = 0; i < 3; i++) {
        m_out << std::setw(5) << grid.npts(i);
        m_out << std::setw(12) << std::setprecision(6) << std::fixed;
        for (size_t j = 0; j < 3; j++) m_out << grid.direction(i, j);
        m_out << std::endl;
    }
    m_out << atoms;
}


} // namespace libwfa



