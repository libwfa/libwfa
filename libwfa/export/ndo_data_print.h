#ifndef LIBWFA_NDO_DATA_PRINT_H
#define LIBWFA_NDO_DATA_PRINT_H

#include <iostream>
#include "ev_printer_i.h"

namespace libwfa {

/** \brief Implementation of ev_printer_i to print NDO summary to output stream

    \ingroup libwfa
 **/
class ndo_data_print : public ev_printer_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    size_t m_nndo; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param nndo # of leading occupation numbers to print
     */
    ndo_data_print(size_t nndo = 3) : m_nndo(nndo) { }

    /** \copydoc ev_printer_i::perform
     **/
    virtual size_t perform(density_type type,
        const ab_vector &ni, std::ostream &out) const;

private:
    size_t print(const arma::Col<double> &ni, std::ostream &out) const;
};


} // namespace libwfa

#endif
