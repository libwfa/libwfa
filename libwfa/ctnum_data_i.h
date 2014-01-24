#ifndef LIBWFA_CTNUM_PRINT_I_H
#define LIBWFA_CTNUM_PRINT_I_H

#include "ab_matrix.h"

namespace libwfa {

/** \brief Interface for printer of CT number data

    \ingroup libwfa
 **/
class ctnum_data_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_data_i() { }

    /** \brief Print the CT number data
        \param om CT number data (omega matrix)
        \param om_tot Omega total
     **/
    virtual void perform(const ab_matrix &om, const double (&om_tot)[2]) = 0;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_PRINT_I_H
