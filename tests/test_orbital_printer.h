#ifndef LIBWFA_TEST_ORBITAL_PRINTER_H
#define LIBWFA_TEST_ORBITAL_PRINTER_H

#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {

/** \brief Orbital printer to test that the right orbitals are printed

    \ingroup libwfa
 **/
class test_orbital_printer : public orbital_printer_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    typedef orbital_type::flag_t ot_flag;

private:
    double m_thresh;
    ot_flag m_ot;

public:
    /** \brief Constructor
        \param out Output stream
        \param title Description of orbitals and densities
        \param dt Flag which density types to export
        \param ot Flag which orbital types to export
     **/
    test_orbital_printer(double thresh = 1e-2,
        const ot_flag &ot = ot_flag(orbital_type::ALL)) :
        m_thresh(thresh), m_ot(ot)
    { }

    /** \brief Destructor
     **/
    virtual ~test_orbital_printer() { }

    /** \copydoc orbital_printer_i::perform

        Appends the orbital type to the header and writes the selected
        orbitals and eigenvalues to the output stream
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s);

    /** \copydoc orbital_printer_i::perform

        Appends the orbital type to the header and writes the selected
        orbitals and eigenvalues to the output stream
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b);

private:
    void check(const arma::vec &e, const orbital_selector &s);
};


} // namespace libwfa

#endif // LIBWFA_TEST_ORBITAL_PRINTER_H
