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
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ctnum_analysis_i &m_analysis; //!< Analysis object
    const ctnum_printer_i &m_printer; //!< Printer

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param a Analysis object
        \param pr Result printer
     **/
    ctnumbers(const arma::Mat<double> &s, const ctnum_analysis_i &a,
        const ctnum_printer_i &pr) : m_s(s), m_analysis(a), m_printer(pr) { }

    /** \brief Destructor
     **/
    ~ctnumbers() { }

    /** \brief Perform analysis
        \param[in] tdm Transition density matrix
        \param[out] om Resulting omega matrix
        \param[out] om_tot Total omega

        Performs the CT number analysis using the analysis object and returns
        the resulting data.
     **/
    void perform(const ab_matrix& tdm,
        ab_matrix &om, double (&om_tot)[2]) const;


    /** \brief Perform analysis
        \param[in] tdm Transition density matrix
        \param[out] out Output stream

        Performs the CT number analysis using the analysis object and prints
        the results using an ctnum_print_i object.
     **/
    void perform(const ab_matrix &tdm, std::ostream &out) const {

        double om_tot[2];
        ab_matrix om;

        perform(tdm, om, om_tot);
        m_printer.perform(om, om_tot, out);
    }
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
