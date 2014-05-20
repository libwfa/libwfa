#ifndef LIBWFA_NO_DATA_PRINT_H
#define LIBWFA_NO_DATA_PRINT_H

#include <iostream>
#include "ev_printer_i.h"

namespace libwfa {

/** \brief Implementation of ev_printer_i to print NO occupation numbers to
        output stream

    \ingroup libwfa
 **/
class no_data_print : public ev_printer_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    std::ostream &m_out; //!< Output stream
    size_t m_nno; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param out Output stream
        \param nno # of leading occupation numbers to print
     */
    no_data_print(std::ostream &out, size_t nno = 3) :
        m_out(out), m_nno(nno) { }

    /** \copydoc ev_printer_i::perform
     **/
    virtual size_t perform(density_type type, const ab_vector &ni);

private:
    size_t print(const arma::Col<double> &ni);
    size_t print_total(const arma::Col<double> &ni);
};


} // namespace libwfa

#endif
