#ifndef LIBWFA_MOLDEN_FILE_BASE_H
#define LIBWFA_MOLDEN_FILE_BASE_H

#include "ab_matrix.h"
#include "ab_vector.h"

namespace libwfa {

/** \brief Molden file object.

    Object representing a molden file.

    TODO: Add additional functions if required.

    \ingroup libwfa
 **/
class molden_file_base {
public:
    virtual ~molden_file_base() { }

    /** \brief Write the orbital coefficients and orbital energies to file.
        \param coeff Orbital coefficients
        \param ene Orbital energies
        \param nocc_a Number of occupied alpha orbitals
        \param nocc_b Number of occupied beta orbitals

        Assume that the first nocc_a (nocc_b) coefficients refer to occupied
        orbitals.
     **/
    virtual void perform(
            const ab_matrix &coeff, const ab_vector &ene,
            size_t nocc_a, size_t nocc_b) = 0;
};


} // namespace adcman

#endif
