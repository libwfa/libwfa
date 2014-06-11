#ifndef LIBWFA_NO_ANALYSIS_H
#define LIBWFA_NO_ANALYSIS_H

#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>

namespace libwfa {


/** \brief Perform complete NO analysis of a state density matrix

     \ingroup libwfa
 **/
class no_analysis {
private:
    const ab_matrix &m_c; //!< MO coefficients
    const ab_matrix &m_sdm; //!< State density matrix
    const ev_printer_i &m_pr; //!< Formating object

public:
    /** \brief Constructor
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param sdm State density matrix
        \param pr Formating object
     **/
    no_analysis(const ab_matrix &c, const ab_matrix &sdm,
        const ev_printer_i &pr) : m_c(c), m_sdm(sdm), m_pr(pr) { }

    /** \brief Perform NTO analysis
        \param pr Printer of NOs
        \param out Output stream
     **/
    void perform(export_data_i &pr, std::ostream &out) const;
};


} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_H
