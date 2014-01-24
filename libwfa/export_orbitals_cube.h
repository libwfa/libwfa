#ifndef LIBWFA_EXPORT_ORBITALS_CUBE_H
#define LIBWFA_EXPORT_ORBITALS_CUBE_H

#include "export_cube_base.h"
#include "export_orbitals_i.h"

namespace libwfa {

/** \brief Class for exporting orbital data to a cube file.

    \ingroup libwfa
 **/
class export_orbitals_cube : public export_orbitals_i {
private:
    export_cube_base &m_core; //!< Export object into cube file

public:
    /** \brief Constructor
        \param core Export object into cube files
     **/
    export_orbitals_cube(export_cube_base &core) : m_core(core) { }

    /** \brief Destructor
     **/
    virtual ~export_orbitals_cube() { }

    /** \copydoc export_orbitals_i::perform
     **/
    virtual void perform(const ab_matrix &coeff,
            const ab_vector &ene, const ab_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_CUBE_H
