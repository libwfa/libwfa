#ifndef LIBWFA_ORBITAL_PRINTER_NIL_H
#define LIBWFA_ORBITAL_PRINTER_NIL_H

#include "orbital_printer_i.h"

namespace libwfa {


/** \brief No export of orbitals

    \ingroup libwfa
 **/
class orbital_printer_nil : public orbital_printer_i {
public:
    /** \brief Destructor
     **/
    virtual ~orbital_printer_nil() { }

    //! \brief Implementation of orbital_printer_i interface
    //@{
    virtual void perform(orbital_type type, const orbital_data &orb,
            const orbital_selector &s) { }

    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b) { }
    //@}
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_NIL_H
