#ifndef LIBWFA_ORBITAL_PRINTER_BASIC_H
#define LIBWFA_ORBITAL_PRINTER_BASIC_H

#include "orbital_printer_i.h"

namespace libwfa {

/** \brief Print density matrices in AO basis

    \ingroup libwfa
 **/
class orbital_printer_basic : public orbital_printer_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef orbital_type::flag_t ot_flag;

private:
    std::ostream &m_out; //!< Output stream
    std::string m_title; //!< Title to put as header
    ot_flag m_ot; //!< Flag which orbital types are printed

public:
    /** \brief Constructor
        \param out Output stream
        \param title Description of orbitals and densities
        \param dt Flag which density types to export
        \param ot Flag which orbital types to export
     **/
    orbital_printer_basic(std::ostream &out, const std::string &title,
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL)) :
        m_out(out), m_title(title), m_ot(ot)
    { }

    /** \brief Destructor
     **/
    virtual ~orbital_printer_basic() { }

    /** \copydoc orbital_printer_i::perform

        Appends the orbital type to the header and writes the selected
        orbitals and eigenvalues to the output stream
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s);

    /** \copydoc orbital_printer_i::perform

        Appends the orbital type to the header and writes the selected
        orbitals and eigenvalues to the output stream
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b);
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_BASIC_H
