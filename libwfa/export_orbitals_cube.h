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
    export_cube_base &m_core; //!< Object to export orbitals into cube file
    std::string m_prefix; //!< Prefix to use for export

public:
    /** \brief Constructor
        \param core Export  into cube files
     **/
    export_orbitals_cube(export_cube_base &core, const std::string &prefix) :
        m_core(core), m_prefix(prefix) { }

    /** \brief Destructor
     **/
    virtual ~export_orbitals_cube() { }

    /** \brief Set the prefix
     **/
    void set_prefix(const std::string &prefix) { m_prefix = prefix; }

    /** \copydoc export_orbitals_i::perform
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_CUBE_H
