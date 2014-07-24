#ifndef LIBWFA_EXPORT_ORBITALS_MOLDEN_H
#define LIBWFA_EXPORT_ORBITALS_MOLDEN_H

#include "export_data_i.h"
#include "export_molden_i.h"

namespace libwfa {

/** \brief Class for exporting orbital data to Molden file.

    Implements the interface export_data_i to export orbital data to a
    Molden file.

    \ingroup libwfa
 **/
class export_orbitals_molden : public export_data_i {
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
    export_orbitals_molden(export_molden_i &core, const std::string &id,
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL));

    /** \brief Destructor
     **/
    virtual ~export_orbitals_molden() { }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(density_type type, const ab_matrix &dm, bool ab_sep = true,
        size_t spin_tr_d = 0) { }

    /** \copydoc export_data_i::perform(libwfa::orbital_type, libwfa::ab_matrix&, libwfa::ab_vector&, libwfa::ab_orbital_selector&)
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_orbital_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_MOLDEN_H
