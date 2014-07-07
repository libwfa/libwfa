#ifndef LIBWFA_TRANSFORMATIONS_H
#define LIBWFA_TRANSFORMATIONS_H

#include "ab_matrix.h"
#include "ab_vector.h"

namespace libwfa {


/** \brief Forms electron and hole density matrices from a transition density
        matrix
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

    \ingroup libwfa
 **/
void form_eh(const arma::Mat<double> &s, const ab_matrix &tdm,
        ab_matrix &de, ab_matrix &dh);

/** \brief Forms omega matrix from a transition density matrix
    \param[in] s Overlap matrix
    \param[in] tdm Transition density matrix
    \param[out] om Omega matrix

    The function implements the transforms of the transition density matrix
    into the Omega matrix, which is further used to compute the charge transfer
    numbers.
    
    The implementation uses the new formula
    \f[
    \Omega_{\mu\nu} =
    0.5 \times \left(
    (\mathbf{D}\mathbf{S})_{\mu\nu} \times (\mathbf{S}\mathbf{D})_{\mu\nu} +
    D_{\mu\nu} \times (\mathbf{S}\mathbf{D}\mathbf{S})_{\mu\nu}
    \right)
    \f]
     from [JCP(2014), DOI: 10.1063/1.4885819] rather than the original
     formula from [JCTC(2012), 8, 2777].

    The output matrices are reshaped and resized as required.

    \ingroup libwfa
 **/
void form_om(const arma::Mat<double> &s, const ab_matrix &tdm,
        ab_matrix &om);

/** \brief Diagonalizes a density matrix in AO basis
    \param[in] s Overlap matrix
    \param[in] c AO2MO coefficient matrix
    \param[in] dm Density matrix
    \param[out] ev Vector of eigenvalues
    \param[out] u Matrix of eigenvectors in the AO basis

    The function transforms the density matrix into an orthogonal basis
    using overlap matrix S and coefficient matrix C
    \f[
    D^{\text{MO}} = (S C)' D (S C)
    \f]
    [this assumes that the column dimension of the coefficient matrix is
    the AO basis]. The resulting density matrix is then diagonalized. The
    eigenvalues are returned unaltered while the eigenvector matrix is
    transformed using the original coefficent matrix to obtain the
    transformation matrix
    \f[
    U = C \tilde{U}
    \quad\text{ with }\quad
    D^{\text{MO}} \tilde{U} = \tilde{U} \Lambda
    \f]
    This transformation matrix can be used to form the density matrix from the
    eigenvalues via
    \f[
    D = U \Lambda U'
    \f]
    and at the same time forms the coefficient matrix for the new orbitals

    \ingroup libwfa
 **/
void diagonalize_dm(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &dm, ab_vector &ev, ab_matrix &u);


/** \brief Constructs attachement and detachment densities
    \param[in] ev Eigenvalues of the difference density matrix
    \param[in] u Eigenvector matrix in the AO basis
    \param[out] da Attachment density matrix
    \param[out] dd Detachment density matrix

    The function sorts the eigenvectors based on the eigenvalues into vectors
    belonging to the attachment and detachment densities, respectively.
    These are then used to back-transform into density matrices. The transform
    is performed as
    \f[
    D = U \Lambda U'
    \f]

    ATTENTION:
    The transformation matrix u is expected to be an eigenvector matrix as
    returned by the function diagonalize_dm.

    \ingroup libwfa
 **/
void form_ad(const ab_vector &ev, const ab_matrix &u,
        ab_matrix &da, ab_matrix &dd);


/** \brief Constructs attachement and detachment densities from a
        difference density matrix
    \param[in] s Overlap matrix
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
void form_ad(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &dm, ab_matrix &da, ab_matrix &dd);


} // namespace libwfa

#endif // LIBWFA_TRANSFORMATIONS_H
