#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NTO_ANALYSIS_H

#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>
#include <libwfa/core/transformations_dm.h>

namespace libwfa {


/** \brief Performs basic NTO analysis using electron and hole density matrices

     \ingroup libwfa
 **/
class nto_analysis_basic {
private:
    ab_vector m_ee, m_eh; //!< Hole and electron eigenvalues
    ab_matrix m_ue, m_uh; //!< Hole and electron eigenvectors

public:
    /** \brief Constructor

        Diagonalize the hole and particle densities and store the results.

        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param edm Electron density matrix
        \param edm Hole density matrix
     **/
    nto_analysis_basic(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm) :
        m_ee(edm.is_alpha_eq_beta()), m_eh(hdm.is_alpha_eq_beta()),
        m_ue(edm.is_alpha_eq_beta()), m_uh(hdm.is_alpha_eq_beta()) {

        diagonalize_dm(s, c, edm, m_ee, m_ue);
        diagonalize_dm(s, c, hdm, m_eh, m_uh);
    }

    /** \brief Return eigenvalues
     *
     * \param elec electron or hole
     */
    ab_vector &get_eigval(bool elec) {
        if (elec) return m_ee;
        else return m_eh;
    }

    /** \brief Return eigenvectors
     *
     * \param elec electron or hole
     */
    ab_matrix &get_eigvect(bool elec) {
        if (elec) return m_ue;
        else return m_uh;
    }

    /** \brief Perform basic NTO analysis

        Eigenvalue print-out and export of MO file.

        \param evpr Eigenvalue printer
        \param opr Printer of NTOs
        \param out Output stream
     **/
    void perform(const ev_printer_i &evpr, export_data_i &opr,
            std::ostream &out) const;

};


/** \brief Perform complete NTO analysis of a transition density matrix

     \ingroup libwfa
 **/
class nto_analysis {
private:
    const arma::mat &m_s; //!< Overlap matrix
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
    nto_analysis(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &tdm, const ev_printer_i &pr) :
        m_s(s), m_c(c), m_tdm(tdm), m_pr(pr) { }

    /** \brief Perform NTO analysis
        \param[out] edm Electron density matrices
        \param[out] hdm Hole density matrices
        \param opr Printer of NTOs
        \param out Output stream
     **/
    void perform(ab_matrix &edm, ab_matrix &hdm, export_data_i &opr,
        std::ostream &out) const;
};

} // namespace libwfa

#endif
