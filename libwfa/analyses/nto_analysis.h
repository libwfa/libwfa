#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NTO_ANALYSIS_H

#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>

namespace libwfa {


/** \brief Performs basic NTO analysis using electron and hole density matrices

     \ingroup libwfa
 **/
class nto_analysis_basic {
private:
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients
    const ab_matrix &m_edm; //!< Electron density matrix
    const ab_matrix &m_hdm; //!< Hole density matrix
    const ev_printer_i &m_pr; //!< Formating object of NTO summary

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param edm Electron density matrix
        \param edm Hole density matrix
        \param pr Formating object
     **/
    nto_analysis_basic(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm, const ev_printer_i &pr) :
        m_s(s), m_c(c), m_edm(edm), m_hdm(hdm), m_pr(pr) { }

    /** \brief Perform basic NTO analysis
        \param pr Printer of NTOs
        \param out Output stream
     **/
    void perform(export_data_i &pr, std::ostream &out) const;
};


/** \brief Perform complete NTO analysis of a transition density matrix

     \ingroup libwfa
 **/
class nto_analysis {
private:
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients
    const ab_matrix &m_tdm; //!< Transition density matrix
    const ev_printer_i &m_pr; //!< Printer of NTO summary

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param tdm Transition density matrices
        \param pr Formating object
     **/
    nto_analysis(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &tdm, const ev_printer_i &pr) :
        m_s(s), m_c(c), m_tdm(tdm), m_pr(pr) { }

    /** \brief Perform NTO analysis
        \param edm Electron density matrices
        \param hdm Hole density matrices
        \param opr Printer of NTOs
        \param out Output stream
     **/
    void perform(ab_matrix &edm, ab_matrix &hdm, export_data_i &opr,
        std::ostream &out) const;
};

} // namespace libwfa

#endif
