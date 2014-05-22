#ifndef LIBWFA_EXPORT_DENSITIES_I_H
#define LIBWFA_EXPORT_DENSITIES_I_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/density_type.h>

namespace libwfa {

/** \brief Interface for exporting density matrices.

    Interface to implement the export of density matrices.

    \ingroup libwfa
 **/
class export_densities_i {
public:
    /** \brief Destructor
     **/
    virtual ~export_densities_i() { }

    /** \brief Export the density matrix
        \param type Type of density matrix
        \param dm Density matrix
     **/
    virtual void perform(density_type type, const ab_matrix &dm) = 0;
};

} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
