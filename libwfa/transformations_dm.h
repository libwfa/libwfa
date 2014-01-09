#ifndef LIBWFA_TRANSFORMATIONS_H
#define LIBWFA_TRANSFORMATIONS_H

#include "ab_matrix.h"
#include "ab_vector.h"

namespace libwfa {


/** \brief Transforms a transition density matrix into electron and hole
        density matrices
    \param[in] s Overlap matrix
    \param[in] tdm Transition density matrix
    \param[out] de Electron density matrix
    \param[out] dh Hole density matrix

    The function implements the transforms of the transition density matrix
    into electron and hole density matrices
    \f[
    E_{\mu\nu} = \sum_{\eta\xi} T_{\mu\eta} S_{\eta\xi} T_{\nu\xi}
    \f]
    and
    \f[
    H_{\mu\nu} = \sum_{\eta\xi} T_{\eta\mu} S_{\eta\xi} T_{\xi\nu}
    \f]

    The output matrices are reshaped and resized as required.

    TODO: check formulas!!! what is the correct density matrix:
        T^{ao} = C' T C ot T^{ao} = S C' T C S?!

    \ingroup libwfa
 **/
void form_eh(const arma::Mat<double> &s, const ab_matrix &tdm,
        ab_matrix &de, ab_matrix &dh);


/** \brief Diagonalizes the given density matrix in AO basis
    \param[in] c AO2MO coefficient matrix
    \param[in] dm Density matrix
    \param[out] ev Vector of eigenvalues
    \param[out] u New coefficient matrix from AO

    The function transforms the density matrix into an orthogonal basis (the MO
    basis), first. Then the density matrix is diagonalized. The eigenvalues are
    returned unaltered while the eigenvector matrix is transformed with the MO
    coefficents to obtain the transformation matrix from AO basis.

    \ingroup libwfa
 **/
void diagonalize_dm(const ab_matrix &c, const ab_matrix &dm,
        ab_vector &ev, ab_matrix &u);


/** \brief Constructs attachement and detachment densities
    \param[in] ev Eigenvalues of the density matrix
    \param[in] u Eigenvectors of the density matrix
    \param[out] da Attachment density matrix
    \param[out] dd Detachment density matrix

    The function sorts the eigenvectors based on the eigenvalues into vectors
    belonging to the attachment and detachment densities, respectively.
    These are then used to back-transform into density matrices.

    \ingroup libwfa
 **/
void form_ad(const ab_vector &ev, const ab_matrix &u,
        ab_matrix &da, ab_matrix &dd);

/** \brief Constructs attachement and detachment densities from a density matrix
    \param[in] c AO2MO coefficient matrix
    \param[in] dm Density matrix
    \param[out] da Attachment density matrix
    \param[out] dd Detachment density matrix

    The function diagonalizes the density matrix first. Then it sorts the
    eigenvectors based on the eigenvalues into vectors belonging to the
    attachment and detachment densities, respectively. These are then used to
    back-transform into density matrices.

    \ingroup libwfa
 **/
void form_ad(const ab_matrix &c, const ab_matrix &dm,
        ab_matrix &da, ab_matrix &dd);

} // namespace adcman

#endif // LIBWFA_TRANSFORMATIONS_H
