#ifndef LIBWFA_DENSITY_PRINTER_BASIC_H
#define LIBWFA_DENSITY_PRINTER_BASIC_H

#include "density_printer_i.h"

namespace libwfa {

/** \brief Print density matrices in AO basis

    \ingroup libwfa
 **/
class density_printer_basic : public density_printer_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef density_type::flag_t dt_flag;

private:
    std::ostream &m_out; //!< Output stream
    std::string m_title; //!< Title to put as header
    dt_flag m_dt; //!< Flag which density types are printed

public:
    /** \brief Constructor
        \param out Output stream
        \param title Description of orbitals and densities
        \param dt Flag which density types to export
     **/
    density_printer_basic(std::ostream &out, const std::string &title,
        const dt_flag &dt = dt_flag(density_type::DT_ALL)) :
        m_out(out), m_title(title), m_dt(dt)
    { }

    /** \brief Destructor
     **/
    virtual ~density_printer_basic() { }

    /** \copydoc density_printer_i::perform

        Appends the density type to the header and writes the density
        matrix to the output stream
     **/
    virtual void perform(density_type type, const ab_matrix &dm);
};


} // namespace libwfa

#endif // LIBWFA_DENSITY_PRINTER_BASIC_H
