#ifndef LIBWFA_NTO_DATA_PRINT_H
#define LIBWFA_NTO_DATA_PRINT_H

#include <iostream>
#include "ev_printer_i.h"

namespace libwfa {

/** \brief Implementation of ev_printer_i to print NTO summary to output stream

    \ingroup libwfa
 **/
class nto_data_print : public ev_printer_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    double m_thresh; //!< Threshold of important NTOs
    size_t m_nnto; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param out Output stream
        \param thresh Threshold for important NTOs
        \param nnto # of leading occupation numbers to print
     */
    nto_data_print(double thresh = 1e-6, size_t nnto = 3) :
        m_thresh(thresh), m_nnto(nnto) { }

    /** \copydoc ev_printer_i::perform
     **/
    virtual size_t perform(density_type type,
            const ab_vector &ni, std::ostream &out) const;

private:
    size_t print(const arma::Col<double> &ni, std::ostream &out) const;
};

} // namespace libwfa

#endif
