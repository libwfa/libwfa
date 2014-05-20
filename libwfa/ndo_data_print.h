#ifndef LIBWFA_NDO_DATA_PRINT_H
#define LIBWFA_NDO_DATA_PRINT_H

#include <iostream>
#include "ev_data_i.h"

namespace libwfa {

/** \brief Implementation of ev_data_i to print NDO occupation numbers to
        output stream

    \ingroup libwfa
 **/
class ndo_data_print : public ev_data_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    std::ostream &m_out; //!< Output stream
    size_t m_nndo; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param out Output stream
        \param nndo # of leading occupation numbers to print
     */
    ndo_data_print(std::ostream &out, size_t nndo = 3) :
        m_out(out), m_nndo(nndo) { }

    /** \copydoc nto_data_i::perform
     **/
    virtual size_t perform(density_type type, const ab_vector &ni);

private:
    size_t print(const arma::Col<double> &ni);
};


} // namespace libwfa

#endif
