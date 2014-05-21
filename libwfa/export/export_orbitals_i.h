#ifndef LIBWFA_EXPORT_ORBITALS_I_H
#define LIBWFA_EXPORT_ORBITALS_I_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/ab_vector.h>
#include <libwfa/core/ab_selector.h>
#include <libwfa/core/orbital_type.h>

namespace libwfa {

/** \brief Interface for exporting orbital data.

    Interface to implement for export of orbital data. The orbital coefficient
    data should be ordered such that the coefficient vectors are the columns
    of the coefficient matrix.

    \ingroup libwfa
 **/
class export_orbitals_i {
public:
    virtual ~export_orbitals_i() { }

    /** \brief Write the given orbital coefficients to the file.
        \param type Orbital type
        \param coeff Orbital coefficients
        \param ev Orbital eigenvalues
        \param s Selector of specific orbital coefficients
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
            const ab_vector &ev, const ab_selector &s) = 0;
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_I_H
