#ifndef LIBWFA_EXPORT_DENSITIES_I_H
#define LIBWFA_EXPORT_DENSITIES_I_H

#include "ab_matrix.h"
#include "density_type.h"

namespace libwfa {

/** \brief No export of density matrices

    \ingroup libwfa
 **/
class export_densities_none : public export_densities_i {
public:
    virtual ~export_densities_none() { }

    /** \copydoc export_densities_i::perform
     **/
    virtual void perform(const std::string &name, const ab_matrix &dm) { }
};

} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_I_H
