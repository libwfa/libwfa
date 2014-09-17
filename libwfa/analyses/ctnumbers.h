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
    const arma::mat &m_s; //!< Overlap matrix
    const ab_matrix &m_tdm; //!< Transition density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param s Overlap matrix
        \param tdm Transition density matrix
     **/
    ctnumbers(const ctnum_analysis_i &a, const arma::mat &s,
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

    /** \brief Forms omega matrix from a transition density matrix
        \param[in] s Overlap matrix
        \param[in] tdm Transition density matrix
        \param[out] om Omega matrix

        The function implements the transforms of the transition density matrix
        into the Omega matrix, which is further used to compute the charge transfer
        numbers.

        The implementation uses the new formula
        \f[
        \Omega_{\mu\nu} = 0.5 \left[
            (\mathbf{D}\mathbf{S})_{\mu\nu} \times (\mathbf{S}\mathbf{D})_{\mu\nu}
            + D_{\mu\nu} \times (\mathbf{S}\mathbf{D}\mathbf{S})_{\mu\nu}
            \right]
        \f]
         from [JCP(2014), DOI: 10.1063/1.4885819] rather than the original
         formula from [JCTC(2012), 8, 2777].

        The output matrices are reshaped and resized as required.
     **/
    static void form_om(const arma::mat &s,
            const ab_matrix &tdm, ab_matrix &om);

private:
    static void form_om(const arma::mat &s,
            const arma::mat &tdm, arma::mat &om) {

        om = 0.5 * ( (tdm * s) % (s * tdm) + tdm % (s * tdm * s) );
    }
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
