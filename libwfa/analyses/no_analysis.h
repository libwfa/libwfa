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
    const ev_printer_i &m_pr; //!< Printer object

public:
    /** \brief Constructor
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param pr Printer object
     **/
    no_analysis(const ab_matrix &c, const ev_printer_i &pr) :
        m_c(c), m_pr(pr) { }

    /** \brief Perform NTO analysis
        \param[in] sdm State density matrix
        \param[out] opr Printer of NOs
        \param[out] out Output stream

        (Assumes the electron density matrix is first in av)
     **/
    void perform(const ab_matrix &sdm,
        export_orbitals_i &opr, std::ostream &out) const;
};

} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_H
