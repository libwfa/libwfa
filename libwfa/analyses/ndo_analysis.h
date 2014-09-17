#ifndef LIBWFA_NDO_ANALYSIS_H
#define LIBWFA_NDO_ANALYSIS_H

#include <iostream>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/orbital_data.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {


/** \brief Performs NDO analysis of a difference density matrix

     \ingroup libwfa
 **/
class ndo_analysis {
private:
    orbital_data *m_ndo[2]; //!< NDO coefficients and occupation numbers

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param ddm Difference density matrix
     **/
    ndo_analysis(const arma::mat &s, const ab_matrix &c, const ab_matrix &ddm);

    /** \brief Destructor
     **/
    ~ndo_analysis();

    /** \brief Forms attachment and detachment densities using the NDOs
        \param at Attachment density matrix
        \param de Detachment density matrix

        The function sorts the NDO eigenvectors based on the eigenvalues
        into vectors belonging to the attachment and detachment densities,
        respectively. These are then used to back-transform into density
        matrices.
     **/
    void form_ad(ab_matrix &at, ab_matrix &de) const;

    /** \brief Perform NDO analysis
        \param out Output stream
        \param nndo Max number of important NDO pairs
        \return Number of important NDO pairs
     **/
    void analyse(std::ostream &out, size_t nndo = 3) const;

    /** \brief Export orbitals
        \param pr Orbital printer
        \param nndo Max number of NDO pairs
     **/
    void export_orbitals(orbital_printer_i &pr, size_t nndo) const;

private:
    static void form_ad(const arma::vec &e, const arma::mat &c,
            arma::mat &at, arma::mat &de);

    static void analysis(std::ostream &out,
            const arma::vec &e, size_t nndo = 3);

    static void bld_selector(const arma::vec &e,
            size_t nndo, orbital_selector &sel);
};

} // namespace libwfa

#endif
