#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NTO_ANALYSIS_H

#include <utility>
#include "ev_data_i.h"
#include "export_densities_i.h"
#include "export_orbitals_i.h"

namespace libwfa {


/** \brief Perform complete NTO analysis of a transition density matrix

     \ingroup libwfa
 **/
class nto_analysis {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param thresh Threshold for NTO selection (default: 1e-6)
        \param nnto Number of leading occ. numbers (default: 3)
     **/
    nto_analysis(const arma::Mat<double> &s, const ab_matrix &c) :
        m_s(s), m_c(c) { }

    /** \brief Perform basic NTO analysis
        \param[in] dm Electron and hole density matrices
        \param[out] u Eigenvectors of electron and hole density
        \param[out] nto_print Printer of NTOs
        \param[out] pr Printer of NTO data

        (Assumes the electron density matrix is first in dm)
     **/
    void perform(const ab_matrix_pair &dm, ab_matrix_pair &u,
            export_orbitals_i &nto_print, ev_data_i &pr) const;

    /** \brief Perform NTO analysis
        \param[in] tdm Transition density matrix
        \param[out] av Average electron and hole density matrices
        \param[out] dm_print Printer of density matrices
        \param[out] nto_print Printer of NTOs
        \param[out] pr Printer of NTO data

        (Assumes the electron density matrix is first in av)
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &dm_print, export_orbitals_i &nto_print,
        ev_data_i &pr) const;
};

} // namespace adcman

#endif
