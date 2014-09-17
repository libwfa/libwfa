#ifndef LIBWFA_EXPORT_DATA_CUBE_H
#define LIBWFA_EXPORT_DATA_CUBE_H

#include "export_cube_i.h"
#include "export_data_i.h"

namespace libwfa {

/** \brief Export density matrices and orbitals as cube files

    Class implementing the interface export_data_i for export of density
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
class export_data_cube : public export_data_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef density_type::flag_t dt_flag;
    typedef orbital_type::flag_t ot_flag;

private:
    export_cube_i &m_core; //!< Core class for export as cube files
    std::string m_id; //!< ID / name of density (or -ies)
    std::string m_desc; //!< Description
    dt_flag m_dt; //!< Flag which density types to export
    ot_flag m_ot; //!< Flag which orbital types to export

public:
    /** \brief Constructor
        \param core Export object into cube files
        \param id ID / name of orbitals and densities
        \param desc Description of orbitals and densities
        \param dt Flag which density types to export
        \param ot Flag which orbital types to export
     **/
    export_data_cube(export_cube_i &core,
        const std::string &id, const std::string &desc,
        const dt_flag &dt = dt_flag(density_type::DT_ALL),
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL)) :
        m_core(core), m_id(id), m_desc(desc), m_dt(dt), m_ot(ot) { }

    /** \brief Destructor
     **/
    virtual ~export_data_cube() { }

    /** \copydoc export_data_i::perform(density_type, const ab_matrix&, bool, size_t)
     **/
    virtual void perform(density_type type, const ab_matrix &dm, bool ab_sep = true,
        size_t spin_tr_d = 0);

    /** \copydoc export_data_i::perform(orbital_type, const ab_matrix&, const ab_vector&, const ab_orbital_selector&)
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_orbital_selector &s);

private:
    void perform(const std::string &name, const std::string &desc,
            const arma::mat &c, const orbital_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DATA_CUBE_H
