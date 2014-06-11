#ifndef LIBWFA_CTNUMBERS_H
#define LIBWFA_CTNUMBERS_H

#include <libwfa/core/ab_matrix.h>
#include "ctnum_analysis_i.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Combines the individual components of the CT number analysis by
    - constructing the \f$\Omega\f$ matrices in AOs from transition density
      matrices
    - computing CT number matrices and the total CT number
      \f$ \Omega_{\text{total}}\f$ .

    \ingroup libwfa
 **/
class ctnumbers {
private:
    const ctnum_analysis_i &m_analysis; //!< Analysis object
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_tdm; //!< Transition density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param s Overlap matrix
        \param tdm Transition density matrix
     **/
    ctnumbers(const ctnum_analysis_i &a, const arma::Mat<double> &s,
        const ab_matrix &tdm) : m_analysis(a), m_s(s), m_tdm(tdm) { }

    /** \brief Destructor
     **/
    ~ctnumbers() { }

    /** \brief Perform analysis
        \param om Resulting omega matrix
        \param om_tot Total omega

        Performs the CT number analysis using the analysis object and returns
        the resulting data.
     **/
    void perform(ab_matrix &om, double (&om_tot)[2]) const;
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
