#ifndef LIBWFA_POP_PRINTER_I_H
#define LIBWFA_POP_PRINTER_I_H

#include "pop_data.h"

namespace libwfa {

/** \brief Interface for printer of population data

    \ingroup libwfa
 **/
class pop_printer_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~pop_printer_i() { }

    /** \brief Print the population data
        \param p Population data
     **/
    virtual void perform(const pop_data &p) = 0;
};

} // namespace libwfa

#endif // LIBWFA_POP_PRINTER_I_H
