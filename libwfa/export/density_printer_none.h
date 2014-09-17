#ifndef LIBWFA_DENSITY_PRINTER_NONE_H
#define LIBWFA_DENSITY_PRINTER_NONE_H

#include "density_printer_i.h"

namespace libwfa {


/** \brief No export of densities

    \ingroup libwfa
 **/
class density_printer_none : public density_printer_i {
public:
    /** \brief Destructor
     **/
    virtual ~density_printer_none() { }

    //! \brief Implementation of density_printer_i interface
    //@{
    virtual void perform(density_type type, const ab_matrix &dm) { }
    //@}
};


} // namespace libwfa

#endif // LIBWFA_DENSITY_PRINTER_NONE_H
