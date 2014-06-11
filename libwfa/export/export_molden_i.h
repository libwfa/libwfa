#ifndef LIBWFA_EXPORT_MOLDEN_I_H
#define LIBWFA_EXPORT_MOLDEN_I_H

#include <string>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/ab_vector.h>

namespace libwfa {

/** \brief Molden file object.

    Object representing a molden file.

    TODO: Add additional functions if required.

    \ingroup libwfa
 **/
class export_molden_i {
public:
    virtual ~export_molden_i() { }

    /** \brief Write the orbital coefficients and orbital energies to file.
        \param name Name associated with the orbitals (use as file prefix!?)
        \param coeff Orbital coefficients
        \param ene Orbital energies
        \param nocc_a Number of occupied alpha orbitals
        \param nocc_b Number of occupied beta orbitals

        Assume that the first nocc_a (nocc_b) coefficients refer to occupied
        orbitals.
     **/
    virtual void perform(const std::string &name, const ab_matrix &coeff,
        const ab_vector &ene, size_t nocc_a, size_t nocc_b) = 0;
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_MOLDEN_I_H
