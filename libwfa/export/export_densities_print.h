#ifndef LIBWFA_EXPORT_DENSITIES_PRINT_H
#define LIBWFA_EXPORT_DENSITIES_PRINT_H

#include "export_densities_i.h"

namespace libwfa {

/** \brief Print density matrices in AO basis

    \ingroup libwfa
 **/
class export_densities_print : public export_densities_i {
private:
    std::ostream &m_out; //!< Output stream
    std::string m_title; //!< Title to put as header

public:
    export_densities_print(std::ostream &out, const std::string &title) :
        m_out(out), m_title(title) { }

    /** \brief Destructor
     **/
    virtual ~export_densities_print() { }

    /** \copydoc export_densities_i::perform

        Appends the density type to the header and writes the density
        matrix to the output stream
     **/
    virtual void perform(density_type type, const ab_matrix &dm);
};

} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_PRINT_H
