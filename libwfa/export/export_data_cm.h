#ifndef LIBWFA_EXPORT_DENSITIES_CM_H
#define LIBWFA_EXPORT_DENSITIES_CM_H

#include "export_data_cube.h"
#include "export_orbitals_molden.h"

namespace libwfa {

/** \brief Export density matrices as cube and orbitals as molden files

    Implementation uses \c export_data_cube and \c export_orbitals_molden

    \sa export_data_cube, export_orbitals_molden

    \ingroup libwfa
 **/
class export_data_cm : public export_data_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef density_type::flag_t dt_flag;
    typedef orbital_type::flag_t ot_flag;

private:
    export_data_cube m_edc; //!< Export as cube
    export_orbitals_molden m_edm; //!< Export as molden

public:
    /** \brief Constructor
        \param ccore Export object into cube files
        \param mcore Export object into molden files
        \param id ID / name of orbitals and densities
        \param desc Description of densities
        \param nbf Total number of basis functions
        \param no_a Number of occupied alpha orbitals
        \param no_b Number of occupied beta orbitals
     **/
    export_data_cm(export_cube_base &ccore, export_molden_i &mcore,
        const std::string &id, const std::string &desc,
        const dt_flag &dt = dt_flag(density_type::DT_ALL),
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL)) :
        m_edc(ccore, id, desc, dt, ot_flag()),
        m_edm(mcore, id, ot) {

    }

    /** \brief Destructor
     **/
    virtual ~export_data_cm() { }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(density_type type, const ab_matrix &dm) {
        m_edc.perform(type, dm);
    }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_orbital_selector &s) {
        m_edm.perform(type, coeff, ev, s);
    }
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DATA_CM_H
