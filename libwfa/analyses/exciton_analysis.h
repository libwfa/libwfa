#ifndef LIBWFA_EXCITON_ANALYSIS_H
#define LIBWFA_EXCITON_ANALYSIS_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>
#include "exciton_moments.h"


namespace libwfa {

using namespace arma;

/** \brief Computes the exciton moments based on a TDM

    Computes the exciton moments using the transition density matrix.

    \ingroup libwfa
 **/
class exciton_analysis {
private:
    const mom_builder_i &m_bld; //!< Builder for multipole moments
    const ab_matrix &m_tdm; //!< Transition density matrix
    size_t m_mmax; //!< Max multipole moment to compute


public:
    exciton_analysis(const mom_builder_i &bld, const ab_matrix &tdm);

    /** \brief Perform analysis
        \param mom Resulting exciton moments
     **/
    void perform(ab_exciton_moments &mom) const;

private:
    void calculate(const arma::mat &tdm, exciton_moments &mom) const;

}; // class exciton_analysis

} // namespace libwfa


#endif
