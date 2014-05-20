#ifndef LIBWFA_EV_PRINTER_I_H
#define LIBWFA_EV_PRINTER_I_H

#include "ab_vector.h"
#include "density_type.h"

namespace libwfa {

/** \brief Interface to export / print eigenvalues of density matrices

    \ingroup libwfa
 **/
class ev_printer_i {
public:
    virtual ~ev_printer_i() { }

    /** \brief Print the eigenvalues
        \param type Type of density matrix
        \param ni Occupation number vector (= vector of eigenvalues)
        \return Number of important eigenvalues
     **/
    virtual size_t perform(density_type type, const ab_vector &ni) = 0;
};


} // namespace libwfa

#endif // LIBWFA_EV_PRINTER_I_H
