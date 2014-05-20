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

    /** \brief Perform CT number analysis for given density
        \param[in] om_ao Omega matrix in AO basis
        \param[out] om_at Resulting CT number matrix per atom

        Routine to perform the CT number analysis for one density matrix. The
        routine should adjust the output matrix to the correct size.
     **/
    virtual void perform(const arma::Mat<double> &om_ao,
            arma::Mat<double> &om_at) const = 0;
};


} // namespace libwfa

#endif
