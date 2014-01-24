#ifndef LIBWFA_CTNUMBERS_H
#define LIBWFA_CTNUMBERS_H

#include "ab_matrix.h"
#include "ctnum_analysis_i.h"
#include "ctnum_print_i.h"
#include "state_info.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Puts together the components of CT number analysis to create
    \f$\omega\f$ matrices and compute \f$\omega_{\text{total}}\f$.

    \ingroup libwfa
 **/
class ctnumbers {
private:
    const ctnum_analysis_i &m_analysis; //!< Class to perform the analysis
    const arma::Mat<double> &m_s; //!< Overlap matrix

public:
    /** \brief Constructor
     **/
    ctnumbers(const ctnum_analysis_i &a, const arma::Mat<double> &s) :
        m_analysis(a), m_s(s) {   }

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
        ab_matrix &om, std::vector<double> &om_tot);


    /** \brief Perform analysis
        \param[in] tdm Transition density matrix
        \param[out] out Output stream
        \param[out] pr Printer for omega data

        Performs the CT number analysis using the analysis object and prints
        the results to output stream and using the ctnum_print_i object.
     **/
    void perform(const ab_matrix &tdm, std::ostream &out, ctnum_print_i &pr);
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
