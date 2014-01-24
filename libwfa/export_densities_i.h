#ifndef LIBWFA_EXPORT_DENSITIES_I_H
#define LIBWFA_EXPORT_DENSITIES_I_H

#include "ab_matrix.h"
#include "dm_type.h"
#include "state_info.h"

namespace libwfa {

/** \brief Interface for exporting density matrices.

    Interface to implement the export of density matrices.

    \ingroup libwfa
 **/
class export_densities_i {
public:
    virtual ~export_densities_i() { }

    /** \brief Export the density matrix
        \param si State information
        \param type Type of density matrix
        \param dm Density matrix
     **/
    virtual void perform(const state_info &si,
            dm_type type, const ab_matrix &dm) = 0;
};

} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
