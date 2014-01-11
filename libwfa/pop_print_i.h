#ifndef LIBWFA_POP_PRINT_I_H
#define LIBWFA_POP_PRINT_I_H

#include "pop_data.h"

namespace libwfa {

/** \brief Interface for printer of population data

    \ingroup libwfa
 **/
class pop_print_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~pop_print_i() { }

    /** \brief Print the population data
        \param p Population data
     **/
    virtual void perform(const pop_data &p) = 0;
};

} // namespace libwfa

#endif // LIBWFA_POP_PRINT_I_H
