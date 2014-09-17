#ifndef LIBWFA_EV_PRINTER_NO_H
#define LIBWFA_EV_PRINTER_NO_H

#include <iostream>
#include "ev_printer_i.h"

namespace libwfa {

/** \brief Implementation of ev_printer_i to print NO summaryy to output stream

    \ingroup libwfa
 **/
class ev_printer_no : public ev_printer_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    size_t m_nno; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param nno # of leading occupation numbers to print
     */
    ev_printer_no(size_t nno = 3) : m_nno(nno) { }

    /** \copydoc ev_printer_i::perform
     **/
    virtual size_t perform(density_type type,
        const ab_vector &ni, std::ostream &out) const;

private:
    size_t print(const arma::vec &ni, std::ostream &out) const;
    size_t print_total(const arma::vec &ni, std::ostream &out) const;
};


} // namespace libwfa

#endif
