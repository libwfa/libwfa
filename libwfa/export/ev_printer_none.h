#ifndef LIBWFA_EV_PRINTER_NONE_H
#define LIBWFA_EV_PRINTER_NONE_H

#include <iostream>
#include "ev_printer_i.h"

namespace libwfa {

/** \brief Implementation of ev_printer_i to print nothing

    \ingroup libwfa
 **/
class ev_printer_none : public ev_printer_i {
public:
    /** \copydoc ev_printer_i::perform
     **/
    virtual size_t perform(density_type type,
            const ab_vector &ni, std::ostream &out) const { return 0; }
};

} // namespace libwfa

#endif
