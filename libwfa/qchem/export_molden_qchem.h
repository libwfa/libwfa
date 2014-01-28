#ifndef LIBWFA_EXPORT_MOLDEN_QCHEM_H
#define LIBWFA_EXPORT_MOLDEN_QCHEM_H

#include "../export_molden_i.h"

namespace libwfa {

/** \brief Export to Molden file object for Q-Chem.

    \ingroup libwfa
 **/
class export_molden_qchem : public export_molden_i {
public:
    /** \brief Constructor
        \param prefix Prefix of molden file
     **/
    export_molden_qchem();

    /** \brief Virtual destructor
     **/
    virtual ~export_molden_qchem() { }

    /** \brief Write the orbital coefficients and orbital energies to file.
        \param coeff Orbital coefficients
        \param ene Orbital energies
        \param nocc_a Number of occupied alpha orbitals
        \param nocc_b Number of occupied beta orbitals

        Assume that the first nocc_a (nocc_b) coefficients refer to occupied
        orbitals.
     **/
    virtual void perform(const std::string &name, const ab_matrix &coeff,
        const ab_vector &ene, size_t nocc_a, size_t nocc_b);
};


} // namespace adcman

#endif // LIBWFA_EXPORT_MOLDEN_QCHEM_H
