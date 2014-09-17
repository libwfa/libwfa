#ifndef LIBWFA_ORBITAL_PRINTER_CUBE_H
#define LIBWFA_ORBITAL_PRINTER_CUBE_H

#include "export_cube_i.h"
#include "orbital_printer_i.h"

namespace libwfa {

/** \brief Export orbitals as cube files

    Class implementing the interface orbital_printer_i for export of density
    matrices as cube files. The real export is done by the core class
    export_cube_base.

    The core export class uses a name and description (both strings) to
    identify densities. This name is constructed by concatenating the ID with
    the density or orbital type
    \code
    m_name + "_" + type.str()
    \endcode
    In case of an unrestricted calculation also the suffixes "_a" and "_b" are
    appended to distinguish alpha and beta spin parts.

    \ingroup libwfa
 **/
class orbital_printer_cube : public orbital_printer_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef orbital_type::flag_t ot_flag;

private:
    export_cube_i &m_core; //!< Core class for export as cube files
    std::string m_id; //!< ID / name of density (or -ies)
    std::string m_desc; //!< Description
    ot_flag m_ot; //!< Flag which orbital types to export

public:
    /** \brief Constructor
        \param core Export object into cube files
        \param id ID / name of orbitals and densities
        \param desc Description of orbitals and densities
        \param dt Flag which density types to export
        \param ot Flag which orbital types to export
     **/
    orbital_printer_cube(export_cube_i &core,
        const std::string &id, const std::string &desc,
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL)) :
        m_core(core), m_id(id), m_desc(desc), m_ot(ot) { }

    /** \brief Destructor
     **/
    virtual ~orbital_printer_cube() { }

    /** \copydoc orbital_printer_i::perform
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s);

    /** \copydoc orbital_printer_i::perform
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b);

private:
    void perform(const std::string &name, const std::string &desc,
            const arma::mat &c, const orbital_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_CUBE_H
