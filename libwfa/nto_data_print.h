#ifndef LIBWFA_NTO_DATA_PRINT_H
#define LIBWFA_NTO_DATA_PRINT_H

#include <iostream>
#include "ev_data_i.h"

namespace libwfa {

/** \brief Implementation of nto_data_i to print to output stream

    \ingroup libwfa
 **/
class nto_data_print : public ev_data_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    std::ostream &m_out; //!< Output stream
    double m_thresh; //!< Threshold of important NTOs
    size_t m_nnto; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param out Output stream
        \param thresh Threshold for important NTOs
        \param nnto # of leading occupation numbers to print
     */
    nto_data_print(std::ostream &out, double thresh = 1e-6, size_t nnto = 3) :
        m_out(out), m_thresh(thresh), m_nnto(nnto) { }

    /** \copydoc nto_data_i::perform
     **/
    virtual size_t perform(density_type type, const ab_vector &ni);

private:
    size_t print(const arma::Col<double> &ni);
};

} // namespace libwfa

#endif
