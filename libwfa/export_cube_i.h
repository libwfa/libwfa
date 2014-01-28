#ifndef LIBWFA_MOLDEN_FILE_I_H
#define LIBWFA_MOLDEN_FILE_I_H

#include <armadillo>
#include <string>
#include <vector>

namespace libwfa {


/** \brief Interface class to export data as cube files.

    Interface class to export certain data (see data_type for possible data) as
    grid data into cube files.

    \ingroup libwfa
 **/
class export_cube_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~export_cube_i() { }

    /** \brief Evaluate an matrix in AO basis on the grid and export as cube
        \param name Name associated with the matrix
        \param mat Matrix data in AO basis
     **/
    virtual void perform(const std::string &name,
        const arma::Mat<double> &mat) = 0;

    /** \brief Evaluate a set of vectors in AO basis on the grid and export as cube
        \param name Name associated with the vectors (used as prefix)
        \param idx Vector of indexes
        \param vecs Set of vectors in AO basis
     **/
    virtual void perform(const std::string &name,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs) = 0;

};

} // namespace adcman

#endif
