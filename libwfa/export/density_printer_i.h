#ifndef LIBWFA_DENSITY_PRINTER_I_H
#define LIBWFA_DENSITY_PRINTER_I_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/density_type.h>

namespace libwfa {

/** \brief Interface for exporting density and orbital data.

    Interface to implement the export of density and orbital data.

    \ingroup libwfa
 **/
class density_printer_i {
public:
    /** \brief Destructor
     **/
    virtual ~density_printer_i() { }

    /** \brief Export the density matrix
        \param type Type of density matrix
        \param dm Density matrix
     **/
    virtual void perform(density_type type, const ab_matrix &dm) = 0;
};

} // namespace libwfa

#endif // LIBWFA_DENSITY_PRINTER_I_H
