#ifndef LIBWFA_CTNUM_PRINTER_I_H
#define LIBWFA_CTNUM_PRINTER_I_H


namespace libwfa {


class ab_matrix;

/** \brief Interface to export / print CT number data

    \ingroup libwfa
 **/
class ctnum_printer_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_printer_i() { }

    /** \brief Print the CT number data
        \param om CT number data (\f$ \Omega \f$ matrix)
        \param om_tot \f$ \Omega_{\text{total}}\f$
        \param out Output stream
     **/
    virtual void perform(const ab_matrix &om,
        const double (&om_tot)[2], std::ostream &out) const = 0;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_PRINTER_I_H
