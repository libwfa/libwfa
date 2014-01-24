#ifndef LIBWFA_CTNUMBERS_H
#define LIBWFA_CTNUMBERS_H

#include "ab_matrix.h"
#include "ctnum_analysis_i.h"
#include "ctnum_print_i.h"
#include "state_info.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Puts together the components of CT number analysis to create
    \f$\omega\f$ matrices and to print them using an implementation of
    ctnum_print_i. In addition \f$\omega_{\text{total}}\f$ is computed.

    \ingroup libwfa
 **/
class ctnumbers {
private:
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ctnum_analysis_i &m_analysis; //!< Class to perform the analysis
    ctnum_print_i &m_printer; //!< Printer of CT number data
    double m_omega[2]; //!< Total \f$\omega\f$

public:
    /** \brief Constructor
     **/
    ctnumbers(const arma::Mat<double> &s, const ctnum_analysis_i &a,
        ctnum_print_i &pr) : m_s(s), m_analysis(a), m_printer(pr) {

        m_omega[0] = m_omega[1] = 0.0;
    }

    /** \brief Destructor
     **/
    ~ctnumbers() { }

    /** \brief
     **/
    double omega(bool alpha) { return m_omega[(alpha ? 0 : 1)]; }

    /** \brief Perform analysis
        \param state State information
        \param tdm Transition density matrix
     **/
    void perform(const state_info &state, const ab_matrix& tdm);
};

} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_H
