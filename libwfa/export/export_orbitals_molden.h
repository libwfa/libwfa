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

private:
    export_molden_i &m_core; //!< Molden file
    std::string m_id; //!< ID / name of orbitals
    size_t m_norbs[4]; //!< Total number of orbitals

public:
    /** \brief Constructor
        \param core Exporter in Molden format
        \param id ID / name of orbitals
        \param no_a Total number of occupied alpha orbitals
        \param nv_a Total number of virtual alpha orbitals
        \param no_b Total number of occupied beta orbitals
        \param nv_b Total number of virtual beta orbitals
     **/
    export_orbitals_molden(export_molden_i &core, const std::string &id,
        size_t no_a, size_t nv_a, size_t no_b, size_t nv_b);

    /** \brief Destructor
     **/
    virtual ~export_orbitals_molden() { }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(density_type type, const ab_matrix &dm) { }

    /** \copydoc export_data_i::perform
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
        const ab_vector &ev, const ab_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_MOLDEN_H
