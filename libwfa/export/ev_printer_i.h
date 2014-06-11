#ifndef LIBWFA_EV_PRINTER_I_H
#define LIBWFA_EV_PRINTER_I_H

#include <libwfa/core/ab_vector.h>
#include <libwfa/core/density_type.h>

namespace libwfa {

/** \brief Interface to export / print eigenvalues of density matrices

    \ingroup libwfa
 **/
class ev_printer_i {
public:
    virtual ~ev_printer_i() { }

    /** \brief Print the eigenvalues
        \param[in] type Type of density matrix
        \param[in] ni Occupation number vector (= vector of eigenvalues)
        \param[out] out Output stream
        \return Number of important eigenvalues
     **/
    virtual size_t perform(density_type type,
        const ab_vector &ni, std::ostream &out) const = 0;
};


} // namespace libwfa

#endif // LIBWFA_EV_PRINTER_I_H
