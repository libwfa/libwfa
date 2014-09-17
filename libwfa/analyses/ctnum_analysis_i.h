#ifndef LIBWFA_CTNUM_ANALYSIS_I_H
#define LIBWFA_CTNUM_ANALYSIS_I_H

#include <armadillo>
#include <cstdlib>

namespace libwfa {

/** \brief Base class for charge transfer (CT) number analysis.

    \ingroup libwfa
 **/
class ctnum_analysis_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_analysis_i() { }
    
    /** \brief Return the dimension of the CT number matrix (e.g. # atoms)
     **/
    virtual size_t size() const = 0;

    /** \brief Perform CT number analysis for given \f$ \Omega \f$ matrix
        \param[in] om_ao \f$ \Omega \f$ matrix in AO basis
        \param[out] om Result CT number matrix

        Perform the CT number analysis for the supplied \f$ \Omega \f$ matrix.
        in AO basis. Any implementation should adjust the output matrix to the
        correct size.
     **/
    virtual void perform(const arma::mat &om_ao,
            arma::mat &om) const = 0;
};


} // namespace libwfa

#endif
