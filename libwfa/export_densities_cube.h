#ifndef LIBWFA_EXPORT_DENSITIES_CUBE_H
#define LIBWFA_EXPORT_DENSITIES_CUBE_H

#include "export_cube_base.h"
#include "export_densities_i.h"

namespace libwfa {

/** \brief Exporting density matrices as cube files

    Class implementing the interface export_densities_i for export of density
    matrices as cube files. The export is done by the core.

    \ingroup libwfa
 **/
class export_densities_cube : public export_densities_i {
private:
    export_cube_base &m_core; //!< Core class for export as cube files

public:
    /** \brief Constructor
        \param core Export object into cube files
     **/
    export_densities_cube(export_cube_base &core) : m_core(core) { }

    /** \brief Destructor
     **/
    virtual ~export_densities_cube() { }

    /** \copydoc export_densities_i::perform
     **/
    virtual void perform(const std::string &name, const ab_matrix &dm);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
