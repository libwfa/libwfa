#ifndef LIBWFA_EV_DATA_I_H
#define LIBWFA_EV_DATA_I_H

#include "ab_vector.h"
#include "density_type.h"

namespace libwfa {

/** \brief Interface class to analyse eigenvalues of density matrices

    \ingroup libwfa
 **/
class ev_data_i {
public:
    virtual ~ev_data_i() { }

    /** \brief Perform the operation
        \param type Type of density matrix
        \param ni Occupation number vector (= vector of eigenvalues)
        \return Number of important occupation numbers
     **/
    virtual size_t perform(density_type type, const ab_vector &ni) = 0;
};


} // namespace adcman

#endif // LIBWFA_EV_DATA_I_H
