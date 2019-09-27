#ifndef LIBWFA_ORBITAL_DATA_H
#define LIBWFA_ORBITAL_DATA_H

#include <armadillo>
#include <vector>
#include <libwfa/soc.h>

namespace libwfa {


/** \brief Storage of oribtal coefficients and occupation numbers (eigenvalues)

    \ingroup libwfa
 **/

class orbital_data {
private:
    arma::vec m_ev; //!< Occupation numbers
    arma::mat m_coeff; //!< Orbital coefficents
    arma::vec m_ene; //!< Orbital "energies" (optional feature)
    std::vector<soc> m_so1e;  //!< Orbital 1e SO integrals (optional feature)
    std::vector<soc> m_somf;  //!< Orbital mean-field SO integrals (optional feature)

public:
    /** \brief Basic constructor
        \param ev Orbital occupation numbers
        \param coeff Orbital coefficients
     **/
    orbital_data(const arma::vec &ev, const arma::mat &coeff) :
        m_ev(ev), m_coeff(coeff) {}

    /** \brief Basic constructor
        \param ev Orbital occupation numbers
        \param coeff Orbital coefficients
        \param coeff Orbital energies
     **/
     orbital_data(const arma::vec &ev, const arma::mat &coeff, const arma::vec &ene) :
        m_ev(ev), m_coeff(coeff), m_ene(ene) {}

    orbital_data(const arma::vec &ev, const arma::mat &coeff, std::vector<soc>& so1e, std::vector<soc>& somf) :
        m_ev(ev), m_coeff(coeff), m_so1e(so1e), m_somf(somf) {}


    /** \brief Constructor
        \param s Overlap matrix
        \param c Coefficient matrix \f$ C \f$
        \param dm Density matrix \f$ D \f$ (in AO basis)

        This constructor assumes that the density matrix \f$ D \f$ has been
        formed using the coefficient matrix \f$ C \f$ via
        \f$ D = C \tilde{D} C' \f$ and that the relation
        \f$ C' S C = \mathbf{I} \f$ holds.

        Then, it transforms the density matrix back to the original basis via
        \f[
        \tilde{D} = C' S D S C
        \f]
        which is then diagonalized. The eigenvalues are stored unaltered as
        occupation numbers, while the eigenvector matrix \f$ \tilde{U} \f$ is
        transformed into the AO basis again to obtain the transformation matrix
        \f[
        U = C \tilde{U}
        \f]
        which is stored as new coefficient matrix.

        The coefficent matrix \f$ U \f$ can be used to form the original
        density matrix from the occupation numbers via
        \f[
        D = U \Lambda U'
        \f]
        Additionally, the occupation numbers can be obtained from the density
        matrix as
        \f[
        \Lambda = U' S D S U
        \f]
     **/
    orbital_data(const arma::mat &s, const arma::mat &c, const arma::mat &dm);

    /** \brief Return Orbital occupation numbers
     **/
    const arma::vec &get_occ() const { return m_ev; }

    /** \brief Return orbital coefficients
     **/
    const arma::mat &get_coeff() const { return m_coeff; }

    /** \brief Tells whether energies were initialized
     **/
    bool got_ene() { return m_ene.size(); }

    /** \brief Tells whether SOC integrals were initialized
     **/
    bool got_soc() { return m_so1e.size(); }

    /** \brief Return orbital energies
     **/
    const arma::vec &get_ene() const { return m_ene; }

    /** \brief Set orbital energies
     **/
    arma::vec &set_ene() { return m_ene; }

    /** \brief Return 1e SOC integrals on NTOs
     **/
    std::vector<soc> &get_so1e() { return m_so1e; }

    /** \brief Return mean-field SOC integrals on NTOs
     **/
    std::vector<soc> &get_somf() { return m_somf; }
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_DATA_H
