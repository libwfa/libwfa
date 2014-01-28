#ifndef LIBWFA_EXPORT_DENSITIES_CUBE_H
#define LIBWFA_EXPORT_DENSITIES_CUBE_H

#include "export_cube_i.h"
#include "export_densities_i.h"

namespace libwfa {

/** \brief Exporting density matrices as cube files

    Class implementing the interface export_densities_i for export of density
    matrices as cube files. The export is done by the core.

    \ingroup libwfa
 **/
class export_densities_cube : public export_densities_i {
private:
    export_cube_i &m_core; //!< Core class for export as cube files
    std::string m_prefix; //!< Prefix for the densities to export

public:
    /** \brief Constructor
        \param core Export object into cube files
     **/
    export_densities_cube(export_cube_i &core, const std::string &prefix) :
        m_core(core), m_prefix(prefix) { }

    /** \brief Destructor
     **/
    virtual ~export_densities_cube() { }

    /** \brief Set the prefix
     **/
    void set_prefix(const std::string &prefix) {  m_prefix = prefix; }

    /** \copydoc export_densities_i::perform
     **/
    virtual void perform(density_type type, const ab_matrix &dm);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
