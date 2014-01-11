#ifndef LIBWFA_CTNUM_PRINT_I_H
#define LIBWFA_CTNUM_PRINT_I_H

#include "ctnum_data.h"

namespace libwfa {

/** \brief Interface for printer of CT number data

    \ingroup libwfa
 **/
class ctnum_print_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_print_i() { }

    /** \brief Print the CT number data
        \param ct CT number data
     **/
    virtual void perform(const ctnum_data &ct) = 0;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_PRINT_I_H
