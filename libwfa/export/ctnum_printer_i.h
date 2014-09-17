#ifndef LIBWFA_CTNUM_PRINTER_I_H
#define LIBWFA_CTNUM_PRINTER_I_H


namespace libwfa {


class ab_matrix;

/** \brief Interface to export CT number data

    \ingroup libwfa
 **/
class ctnum_printer_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_printer_i() { }

    /** \brief Print the CT number data
        \param out Output stream
        \param om CT number data (\f$ \Omega \f$ matrix)
     **/
    virtual void perform(const ab_matrix &om) const = 0;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_PRINTER_I_H
