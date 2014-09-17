#ifndef LIBWFA_NDO_ANALYSIS_H
#define LIBWFA_NDO_ANALYSIS_H

#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>

namespace libwfa {


/** \brief Performs complete NDO analysis of a difference density matrix

     \ingroup libwfa
 **/
class ndo_analysis {
private:
    const arma::mat &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients
    const ab_matrix &m_ddm; //!< Difference density matrix
    const ev_printer_i &m_pr; //!< Formating object
public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param ddm Difference density matrix
        \param pr Formating object
     **/
    ndo_analysis(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &ddm, const ev_printer_i &pr) :
        m_s(s), m_c(c), m_ddm(ddm), m_pr(pr) { }

    /** \brief Perform NDO analysis
        \param[out] at Attachment density matrix
        \param[out] de Detachment density matrix
        \param[out] u  NDO coefficients
        \param[out] ev NDO eigenvalues
        \param pr Printer of NTOs
        \param out Output stream
     **/
    void perform(ab_matrix &at, ab_matrix &de,
        ab_matrix &u, ab_vector &ev,
        export_data_i &opr, std::ostream &out) const;

    /** \brief Perform NDO analysis
        \param[out] at Attachment density matrix
        \param[out] de Detachment density matrix
        \param pr Printer of NTOs
        \param out Output stream
     **/
    void perform(ab_matrix &at, ab_matrix &de,
        export_data_i &opr, std::ostream &out) const;
        
        
    /** \brief Perform NDO analysis
        \param pr Printer of density matrices and NTOs
        \param out Output stream
     **/
    void perform(export_data_i &opr, std::ostream &out) const;
};

} // namespace libwfa

#endif
