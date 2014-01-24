#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NO_ANALYSIS_H

#include <utility>
#include "ab_matrix.h"
#include "ev_data_i.h"
#include "export_orbitals_i.h"

namespace libwfa {


/** \brief Perform complete NO analysis of a state density matrix

     \ingroup libwfa
 **/
class no_analysis {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const ab_matrix &m_c; //!< MO coefficients

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param thresh Threshold for NTO selection (default: 1e-6)
        \param nnto Number of leading occ. numbers (default: 3)
     **/
    no_analysis(const ab_matrix &c) : m_c(c) { }

    /** \brief Perform NTO analysis
        \param[in] sdm State density matrix
        \param[out] no_print Printer of NOs
        \param[out] pr Printer of NO data

        (Assumes the electron density matrix is first in av)
     **/
    void perform(const ab_matrix &sdm,
        export_orbitals_i &nto_print, ev_data_i &pr) const;
};

} // namespace adcman

#endif
