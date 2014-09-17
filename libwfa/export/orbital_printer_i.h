#ifndef LIBWFA_ORBITAL_PRINTER_I_H
#define LIBWFA_ORBITAL_PRINTER_I_H

#include <libwfa/core/orbital_data.h>
#include <libwfa/core/orbital_selector.h>
#include <libwfa/core/orbital_type.h>

namespace libwfa {

/** \brief Interface for printing orbital data.

    \ingroup libwfa
 **/
class orbital_printer_i {
public:
    /** \brief Destructor
     **/
    virtual ~orbital_printer_i() { }

    /** \brief Write the given orbital coefficients and occupation numbers.
        \param type Orbital type
        \param orb Orbital data
        \param s Selector of specific orbitals

        The orbital coefficient vectors can be expected as columns of the
        coefficient matrix.
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s) = 0;

    /** \brief Write the given orbital coefficients and occupation numbers.
        \param type Orbital type
        \param orb_a \f$ \alpha \f$ orbital data
        \param s_a Selector of specific \f$ \alpha \f$ orbitals
        \param orb_b \f$ \beta \f$ orbital data
        \param s_b Selector of specific \f$ \beta \f$ orbitals

        The orbital coefficient vectors can be expected as columns of the
        coefficient matrix.
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b) = 0;
};

} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_I_H
