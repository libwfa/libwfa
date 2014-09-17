#ifndef LIBWFA_NTO_ANALYSIS_H
#define LIBWFA_NTO_ANALYSIS_H

#include <iostream>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/orbital_data.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {


/** \brief Perform complete NTO analysis of a transition density matrix


    \ingroup libwfa
 **/
class nto_analysis {
private:
    orbital_data *m_nto[4]; //!< NTOs

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param tdm Transition density matrix
     **/
    nto_analysis(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &tdm);

    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param edm Electron density matrix
        \param hdm Hole density matrix
     **/
    nto_analysis(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm);

    /** \brief Is alpha == beta
     **/
    bool is_alpha_eq_beta() const { return m_nto[2] == 0; }

    /** \brief Retrieve NTOs
        \param el True: electron NTOs; false: hole NTOs
        \param spin True: beta NTOs; false: alpha NTOs

        TODO: Is spin definition best this way?
     **/
    const orbital_data &get_ntos(bool el, bool spin) const {
        return *m_nto[(el ? 0 : 1) + (spin && m_nto[2] ? 2 : 0)];
    }

    /** \brief Perform NTO analysis
        \param out Output stream
     **/
    void analyse(std::ostream &out, size_t nnto = 3) const;

    /** \brief Export NTOs
        \param pr Orbital printer
        \param thresh Zero threshold (default: \f$ 10^-6 \f$)

     **/
    void export_orbitals(orbital_printer_i &pr, double thresh = 1e-6) const;

    /** \brief Construct the electron and hole density matrices
        \param s Overlap matrix
        \param tdm Transition density matrix
        \param edm Electron density matrix
        \param hdm Hole density matrix

        The function implements the transforms of the transition density
        matrix into electron and hole density matrices
        \f[
            E_{\mu\nu} = \sum_{\eta\xi} T_{\mu\eta} S_{\eta\xi} T_{\nu\xi}
        \f]
        and
        \f[
            H_{\mu\nu} = \sum_{\eta\xi} T_{\eta\mu} S_{\eta\xi} T_{\xi\nu}
        \f]
        The output matrices are reshaped and resized as required.
     **/
    static void form_eh(const arma::mat &s, const ab_matrix &tdm,
        ab_matrix &edm,  ab_matrix &hdm);

private:
    void initialize(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm);

    static void analysis(std::ostream &out,
        const arma::vec &e, size_t nnto = 3);

    static void build_selector(const arma::vec &e, const arma::vec &h,
        double thresh, orbital_selector &sel);

    static void form_eh(const arma::mat &s, const arma::mat &tdm,
        arma::mat &edm, arma::mat &hdm) {
        edm = tdm.t() * s * tdm;
        hdm = tdm * s * tdm.t();
    }
};

} // namespace libwfa

#endif
