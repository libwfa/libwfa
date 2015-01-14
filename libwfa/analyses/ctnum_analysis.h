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
    const arma::mat &m_s; //!< Overlap matrix
    const arma::uvec &m_b2p; //!< Map of basis functions to fragments

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param b2p Map of basis functions to molecular parts or fragments
     **/
    ctnum_analysis(const arma::mat &s, const arma::uvec &b2p);

    /** \brief Virtual destructor
     **/
    virtual ~ctnum_analysis() { }

    /** \copydoc ctnum_analysis_i::size
     **/
    virtual size_t size() const { return m_nparts; }

    /** \copydoc ctnum_analysis_i::perform
     **/
    virtual void perform(const arma::mat& tdm, arma::mat &om) const;

    /** \brief Forms omega matrix from a transition density matrix
        \param[in] s Overlap matrix
        \param[in] tdm Transition density matrix
        \param[out] om Omega matrix

        The function implements the transforms of the transition density
        matrix into the \f$\Omega\f$ matrix, which is further used to compute
        the charge transfer numbers.

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
            const arma::mat &tdm, arma::mat &om);
};

} // namespace libwfa

#endif
