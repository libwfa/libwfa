#ifndef LIBWFA_NO_ANALYSIS_H
#define LIBWFA_NO_ANALYSIS_H

#include <utility>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_orbitals_i.h>

namespace libwfa {


/** \brief Perform complete NO analysis of a state density matrix

     \ingroup libwfa
 **/
class no_analysis {
private:
    const ab_matrix &m_c; //!< MO coefficients

public:
    /** \brief Constructor
        \param c Orbital coefficient matrix for transform in orthogonal basis
     **/
    no_analysis(const ab_matrix &c) : m_c(c) { }

    /** \brief Perform NTO analysis
        \param[in] sdm State density matrix
        \param[out] pr_o Printer of NOs
        \param[out] pr_e Printer of NO occupation numbers

        (Assumes the electron density matrix is first in av)
     **/
    void perform(const ab_matrix &sdm, export_orbitals_i &pr_o,
        ev_printer_i &pr_e) const;
};

} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_H
