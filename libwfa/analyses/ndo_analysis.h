#ifndef LIBWFA_NDO_ANALYSIS_H
#define LIBWFA_NDO_ANALYSIS_H

#include <utility>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_densities_i.h>
#include <libwfa/export/export_orbitals_i.h>

namespace libwfa {


/** \brief Perform complete NDO analysis of a transition density matrix

     \ingroup libwfa
 **/
class ndo_analysis {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients
    const ev_printer_i &m_pr; //!< Printer object
public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param pr Printer object
     **/
    ndo_analysis(const arma::Mat<double> &s, const ab_matrix &c,
        ev_printer_i &pr) : m_s(s), m_c(c), m_pr(pr) { }

    /** \brief Perform NDO analysis
        \param[in] ddm Difference density matrix
        \param[out] ad Attachment / detachment densities (first attach)
        \param[out] opr Printer of NTOs
        \param[out] out Printer object
     **/
    void perform(const ab_matrix &ddm, ab_matrix_pair &ad,
        export_orbitals_i &opr, std::ostream &out) const;

    /** \brief Perform NDO analysis
        \param[in] ddm Difference density matrix
        \param[out] dpr Printer of density matrices
        \param[out] opr Printer of NTOs
        \param[out] out Printer object
     **/
    void perform(const ab_matrix &ddm, export_densities_i &dpr,
        export_orbitals_i &opr, std::ostream &out) const;
};

} // namespace libwfa

#endif
