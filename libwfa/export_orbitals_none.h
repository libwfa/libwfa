#ifndef LIBWFA_EXPORT_ORBITALS_NONE_H
#define LIBWFA_EXPORT_ORBITALS_NONE_H

#include "ab_matrix.h"
#include "density_type.h"

namespace libwfa {

/** \brief No export of orbitals

    \ingroup libwfa
 **/
class export_orbitals_none : public export_orbitals_i {
public:
    virtual ~export_orbitals_none() { }

    /** \copydoc export_orbitals_i::perform
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
            const ab_vector &ev, const ab_selector &s) { }
};

} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_NONE_H
