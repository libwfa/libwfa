#ifndef LIBWFA_DENSITY_PRINTER_NIL_H
#define LIBWFA_DENSITY_PRINTER_NIL_H

#include "density_printer_i.h"

namespace libwfa {


/** \brief No export of densities

    \ingroup libwfa
 **/
class density_printer_nil : public density_printer_i {
public:
    /** \brief Destructor
     **/
    virtual ~density_printer_nil() { }

    //! \brief Implementation of density_printer_i interface
    //@{
    virtual void perform(density_type type, const ab_matrix &dm) { }
    //@}
};


} // namespace libwfa

#endif // LIBWFA_DENSITY_PRINTER_NIL_H
