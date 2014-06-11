#ifndef LIBWFA_EXPORT_DATA_PRINT_H
#define LIBWFA_EXPORT_DATA_PRINT_H

#include "export_data_i.h"

namespace libwfa {

/** \brief Print density matrices in AO basis

    \ingroup libwfa
 **/
class export_data_print : public export_data_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    std::ostream &m_out; //!< Output stream
    std::string m_title; //!< Title to put as header

public:
    export_data_print(std::ostream &out, const std::string &title) :
        m_out(out), m_title(title) { }

    /** \brief Destructor
     **/
    virtual ~export_data_print() { }

    /** \copydoc export_data_i::perform

        Appends the density type to the header and writes the density
        matrix to the output stream
     **/
    virtual void perform(density_type type, const ab_matrix &dm);

    /** \copydoc export_data_i::perform

        Appends the orbital type to the header and writes the selected
        orbitals and eigenvalues to the output stream
     **/
    virtual void perform(orbital_type type, const ab_matrix &coeff,
            const ab_vector &ev, const ab_selector &s);

};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DATA_PRINT_H
