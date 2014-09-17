#ifndef LIBWFA_ORBITAL_PRINTER_MOLDEN_H
#define LIBWFA_ORBITAL_PRINTER_MOLDEN_H

#include "orbital_printer_i.h"
#include "export_molden_i.h"

namespace libwfa {

/** \brief Class for exporting orbital data to Molden file.

    Implements the interface export_data_i to export orbital data to a
    Molden file.

    \ingroup libwfa
 **/
class orbital_printer_molden : public orbital_printer_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef orbital_type::flag_t ot_flag;

private:
    export_molden_i &m_core; //!< Molden file
    std::string m_id; //!< ID / name of orbitals
    ot_flag m_ot; //!< Flag which orbital types are printed

public:
    /** \brief Constructor
        \param core Exporter in Molden format
        \param id ID / name of orbitals
        \param ot Flag which orbital types to export
     **/
    orbital_printer_molden(export_molden_i &core, const std::string &id,
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL));

    /** \brief Destructor
     **/
    virtual ~orbital_printer_molden() { }

    /** \copydoc orbital_printer_i::perform
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s);

    /** \copydoc orbital_printer_i::perform
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b);
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_MOLDEN_H
