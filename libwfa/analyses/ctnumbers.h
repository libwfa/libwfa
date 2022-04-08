#ifndef LIBWFA_CTNUMBERS_H
#define LIBWFA_CTNUMBERS_H

#include <iostream>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/export/ctnum_printer_i.h>
#include "ctnum_analysis_i.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Combines the individual components of the CT number analysis by
    - constructing the \f$\Omega\f$ matrices in AOs from transition density
      matrices
    - computing CT number matrices and the total CT number
      \f$ \Omega_{\mathrm{total}}\f$ .

    \ingroup libwfa
 **/
class ctnumbers {
private:
    double m_tot[2]; //!< Total omega
    double m_om_ab;  //!< Frobenius scalar product of alpha and beta 1TDM
    double m_Phe[2]; //!< Electron-hole permutation
    ab_matrix m_om; //!< Omega matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param tdm Transition density matrix
     **/
    ctnumbers(const ctnum_analysis_i &a, const ab_matrix &tdm);

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

private:
    void initialize(const ctnum_analysis_i &a, const ab_matrix &om);
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
