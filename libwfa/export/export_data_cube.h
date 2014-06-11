#ifndef LIBWFA_EXPORT_DENSITIES_CUBE_H
#define LIBWFA_EXPORT_DENSITIES_CUBE_H

#include "export_cube_base.h"
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

private:
    export_cube_base &m_core; //!< Core class for export as cube files
    std::string m_id; //!< ID / name of density (or -ies)
    std::string m_desc; //!< Description

public:
    /** \brief Constructor
        \param core Export object into cube files
        \param id ID / name of orbitals and densities
        \param desc Description of orbitals and densities
     **/
    export_data_cube(export_cube_base &core, const std::string &id,
        const std::string &desc) :
        m_core(core), m_id(id), m_desc(desc) { }

    /** \brief Destructor
     **/
    virtual ~export_data_cube() { }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(density_type type, const ab_matrix &dm);

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_selector &s);

private:
    void perform(const std::string &name, const std::string &desc,
            const arma::Mat<double> &c, const selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DATA_CUBE_H
