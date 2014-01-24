#ifndef LIBWFA_EXPORT_ORBITALS_MOLDEN_H
#define LIBWFA_EXPORT_ORBITALS_MOLDEN_H

#include "export_orbitals_i.h"
#include "molden_file_base.h"

namespace libwfa {

/** \brief Class for exporting orbital data to Molden file.

    Implements the interface export_orbitals_i to export orbital data to a
    Molden file.

    \ingroup libwfa
 **/
class export_orbitals_molden : public export_orbitals_i {
private:
    molden_file_base &m_file; //!< Molden file
    size_t m_norbs[4]; //!< Total number of orbitals

public:
    /** \brief Constructor
        \param file Molden file
        \param no_a Total number of occupied alpha orbitals
        \param nv_a Total number of virtual alpha orbitals
        \param no_b Total number of occupied beta orbitals
        \param nv_b Total number of virtual beta orbitals
     **/
    export_orbitals_molden(molden_file_base &file,
        size_t no_a, size_t nv_a, size_t no_b, size_t nv_b);

    /** \brief Destructor
     **/
    virtual ~export_orbitals_molden() { }

    /** \copydoc export_orbitals_i::perform
     **/
    virtual void perform(const ab_matrix &coeff,
            const ab_vector &ene, const ab_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_MOLDEN_H
