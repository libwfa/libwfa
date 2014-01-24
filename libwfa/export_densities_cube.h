#ifndef LIBWFA_EXPORT_DENSITIES_CUBE_H
#define LIBWFA_EXPORT_DENSITIES_CUBE_H

#include "export_cube_base.h"
#include "export_densities_i.h"

namespace libwfa {

/** \brief Interface for exporting density matrices.

    Interface to implement the export of density matrices.

    \ingroup libwfa
 **/
class export_densities_cube : public export_densities_i {
private:
    export_cube_base &m_core; //!< Export class into cube files

public:
    export_densities_cube(export_cube_base &core) : m_core(core) { }

    virtual ~export_densities_cube() { }

    /** \copydoc export_densities_i::perform(dm_type, size_t, const ab_matrix &)
     **/
    virtual void perform(dm_type type, size_t idx, const ab_matrix &dm);

    /** \copydoc export_densities_i::perform(const dm_list &)
     **/
    virtual void perform(const dm_list &lst);

private:
    static cube::data_type determine_data_type(dm_type type, bool alpha);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
