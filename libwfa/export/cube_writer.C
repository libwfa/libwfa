#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "cube_writer.h"

namespace libwfa {


cube_writer::cube_writer(const std::string &filename,
    const std::string &line1, const std::string &line2,
    const grid3d &grid, const arma::uvec &atnum, const arma::mat &coords) :
    m_out(filename.c_str()), m_npoints(grid.size()),
    m_nrec(grid.npts().at(2)), m_cur(0) {

    grid.check();
    write_header(line1, line2, grid, atnum, coords);
}


size_t cube_writer::write(const arma::vec &data) {

    // Set the correct output format for data
    m_out << std::setprecision(5) << std::scientific;

    const double *ptr = data.mem;
    size_t i = 0, j = m_cur % m_nrec, k = j % 6;
    while (i < data.n_elem && m_cur < m_npoints) {
        m_out << " " << std::setw(13) << *ptr++;
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


void cube_writer::write_header(const std::string &l1, const std::string &l2,
    const grid3d &grid, const arma::uvec &atnum, const arma::mat &coords) {

    if (atnum.n_elem != coords.n_cols) {
        throw libwfa_exception("cube_writer", "write_header(...)",
                __FILE__, __LINE__, "# atoms");
    }

    m_out << l1 << std::endl;
    m_out << l2 << std::endl;
    m_out << std::setw(5) << atnum.n_elem;

    m_out << std::setprecision(6) << std::fixed;

    arma::subview<double> x0 = grid.origin();
    for (size_t i = 0; i < 3; i++) m_out << std::setw(12) << x0(i);
    m_out << std::endl;

    const arma::Col<unsigned int> &np = grid.npts();
    for (size_t i = 0; i < 3; i++) {
        m_out << std::setw(5) << np(i);

        arma::subview<double> ei = grid.direction(i);
        m_out << std::setprecision(6) << std::fixed;
        for (size_t j = 0; j < 3; j++) m_out << std::setw(12) << ei(j);
        m_out << std::endl;
    }

    for (size_t i = 0; i < atnum.n_elem; i++) {
        m_out << std::setw(5) << atnum(i);
        m_out << std::setprecision(6) << std::fixed;
        m_out << std::setw(12) << 1.0 * atnum(i);
        for (size_t j = 0; j < 3; j++)
            m_out << std::setw(12) << coords(j, i);
        m_out << std::endl;
    }
}


} // namespace libwfa



