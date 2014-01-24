#ifndef LIBWFA_EXPORT_DENSITIES_I_H
#define LIBWFA_EXPORT_DENSITIES_I_H

#include "ab_matrix.h"
#include "dm_list.h"

namespace libwfa {

/** \brief Interface for exporting density matrices.

    Interface to implement the export of density matrices.

    \ingroup libwfa
 **/
class export_densities_i {
public:
    virtual ~export_densities_i() { }

    /** \brief Export the density matrix
        \param type Type of density matrix
        \param idx Index of density matrix
        \param dm Density matrix
     **/
    virtual void perform(dm_type type, size_t idx, const ab_matrix &dm) = 0;

    /** \brief Export the list of density matrices
        \param lst Density matrix list
     **/
    virtual void perform(const dm_list &lst) = 0;
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
