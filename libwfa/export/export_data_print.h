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

public:
    typedef density_type::flag_t dt_flag;
    typedef orbital_type::flag_t ot_flag;

private:
    std::ostream &m_out; //!< Output stream
    std::string m_title; //!< Title to put as header
    dt_flag m_dt; //!< Flag which density types are printed
    ot_flag m_ot; //!< Flag which orbital types are printed

public:
    /** \brief Constructor
        \param out Output stream
        \param title Description of orbitals and densities
        \param dt Flag which density types to export
        \param ot Flag which orbital types to export
     **/
    export_data_print(std::ostream &out, const std::string &title,
        const dt_flag &dt = dt_flag(density_type::DT_ALL),
        const ot_flag &ot = ot_flag(orbital_type::OT_ALL)) :
        m_out(out), m_title(title), m_dt(dt), m_ot(ot)
    { }

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
