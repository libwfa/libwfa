#ifndef LIBWFA_EXPORT_DATA_NONE_H
#define LIBWFA_EXPORT_DATA_NONE_H

#include "export_data_i.h

namespace libwfa {


/** \brief No export of orbitals and densities

    \ingroup libwfa
 **/
class export_data_none : public export_data_i {
public:
    /** \brief Destructor
     **/
    virtual ~export_densities_none() { }

    //! \brief Implementation of export_data_i interface
    //@{
    virtual void perform(density_type type, const ab_matrix &dm) { }

    virtual void perform(orbital_type type, const ab_matrix &coeff,
            const ab_vector &ev, const ab_selector &s) { }
//@}
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DATA_NONE_H
