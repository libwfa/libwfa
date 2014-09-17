#ifndef LIBWFA_CTNUMBERS_H
#define LIBWFA_CTNUMBERS_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/export/ctnum_printer_i.h>
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
    double m_tot[2]; //!< Total omega
    ab_matrix m_om; //!< Omega matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param s Overlap matrix
        \param tdm Transition density matrix
     **/
    ctnumbers(const ctnum_analysis_i &a,
        const arma::mat &s, const ab_matrix &tdm);

    /** \brief Constructor
        \param a Analysis object
        \param om Omega matrix
     **/
    ctnumbers(const ctnum_analysis_i &a, const ab_matrix &om) :
        m_om(om.is_alpha_eq_beta()) {
        initialize(a, om);
    }

    /** \brief Destructor
     **/
    ~ctnumbers() { }

    /** \brief Return omega total
        \param spin If true: beta spin, false: alpha spin
     **/
    const double &omega_total(bool spin) const {
        return m_tot[(spin && ! m_om.is_alpha_eq_beta() ? 1 : 0)];
    }

    /** \brief Return omega
     **/
    const ab_matrix &omega() const { return m_om; }

    /** \brief Perform analysis
        \param out Output stream
     **/
    void analyse(std::ostream &out) const;

    /** \brief Export CT numbers
        \param pr CT numbers printer
     **/
    void do_export(ctnum_printer_i &pr) const { pr.perform(m_om); }

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
            const ab_matrix &tdm, ab_matrix &om);

private:
    void initialize(const ctnum_analysis_i &a, const ab_matrix &om);

    static void form_om(const arma::mat &s,
            const arma::mat &tdm, arma::mat &om) {

        om = 0.5 * ( (tdm * s) % (s * tdm) + tdm % (s * tdm * s) );
    }

};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
