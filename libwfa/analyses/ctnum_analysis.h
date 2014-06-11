#ifndef LIBWFA_CTNUM_ANALYSIS_H
#define LIBWFA_CTNUM_ANALYSIS_H

#include <vector>
#include "ctnum_analysis_i.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Basic implementation of CT number analysis. Provided a mapping of basis
    functions to molecular fragments (e.g. atoms) the class computes for every
    \f$ \Omega \f$ matrix in AOs the CT number matrix in terms of the fragments.

    \ingroup libwfa
 **/
class ctnum_analysis : public ctnum_analysis_i {
private:
    size_t m_nparts; //!< Number of parts
    const arma::Col<size_t> &m_b2p; //!< Map of basis functions to fragments

public:
    /** \brief Constructor
        \param b2c Map of basis functions to molecular parts or fragments
     **/
    ctnum_analysis(const arma::Col<size_t> &b2p);

    virtual ~ctnum_analysis() { }

    /** \copydoc ctnum_analysis_i::size
     **/
    virtual size_t size() const { return m_nparts; }

    /** \copydoc ctnum_analysis_i::perform
     **/
    virtual void perform(const arma::Mat<double>& om_ao,
            arma::Mat<double> &om) const;
};

} // namespace libwfa

#endif
