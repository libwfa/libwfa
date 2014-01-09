#ifndef LIBWFA_POP_PRINT_I_H
#define LIBWFA_POP_PRINT_I_H

#include "pop_data.h"

namespace libwfa {


/** \brief Interface for printer of population data

    \ingroup libwfa
 **/
class pop_print_i {
public:
    virtual ~pop_print_i() { }

    virtual void perform(const pop_data &p) = 0;
};


} // namespace adcman

#endif // LIBWFA_POP_PRINT_I_H
