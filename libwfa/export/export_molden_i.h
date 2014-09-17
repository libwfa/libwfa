#ifndef LIBWFA_EXPORT_MOLDEN_I_H
#define LIBWFA_EXPORT_MOLDEN_I_H

#include <string>
#include <armadillo>

namespace libwfa {

/** \brief Molden file object.

    Object representing a molden file.

    TODO: Add additional functions if required.

    \ingroup libwfa
 **/
class export_molden_i {
public:
    virtual ~export_molden_i() { }

    /** \brief Write the orbital coefficients and orbital energies to file
            (restricted version)
        \param name Name associated with the orbitals (use as file prefix!?)
        \param c Orbital coefficients
        \param e Orbital energies
        \param nocc Number of occupied orbitals

        Assume that the first nocc coefficients refer to occupied orbitals.
     **/
    virtual void perform(const std::string &name, const arma::mat &c,
        const arma::vec &e, size_t nocc) = 0;

    /** \brief Write the orbital coefficients and orbital energies to file.
            (unrestricted version)
        \param name Name associated with the orbitals (use as file prefix!?)
        \param c_a Alpha orbital coefficients
        \param e_a Alpha orbital energies
        \param nocc_a Number of occupied alpha orbitals
        \param c_b Beta orbital coefficients
        \param e_b Beta orbital energies
        \param nocc_b Number of occupied beta orbitals

        Assume that the first nocc_a (nocc_b) coefficients refer to occupied
        orbitals.
     **/
    virtual void perform(const std::string &name,
        const arma::mat &c_a, const arma::vec &e_a, size_t nocc_a,
        const arma::mat &c_b, const arma::vec &e_b, size_t nocc_b) = 0;
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_MOLDEN_I_H
