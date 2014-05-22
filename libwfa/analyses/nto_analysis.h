#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NTO_ANALYSIS_H

#include <utility>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_orbitals_i.h>

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
    const ev_printer_i &m_pr; //!< Printer of NTO summary

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
     **/
    nto_analysis(const arma::Mat<double> &s, const ab_matrix &c,
        const ev_printer_i &pr) : m_s(s), m_c(c), m_pr(pr) { }

    /** \brief Perform basic NTO analysis
        \param[in] eh Particle and hole density matrices (particle is first)
        \param[out] u Eigenvectors of electron and hole density
        \param[out] pr_o Printer of NTOs
        \param[out] out Output stream

        (Assumes the electron density matrix is first in dm)
     **/
    void perform(const ab_matrix_pair &eh, ab_matrix_pair &u,
            export_orbitals_i &pr_o, std::ostream &out) const;

    /** \brief Perform NTO analysis
        \param[in] tdm Transition density matrix
        \param[out] eh Particle and hole density matrices (particle is first)
        \param[out] pr_o Printer of NTOs
        \param[out] out Output stream
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &eh,
        export_orbitals_i &pr_o, std::ostream &out) const;
};

} // namespace libwfa

#endif
