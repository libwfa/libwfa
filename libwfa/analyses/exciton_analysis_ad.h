#ifndef LIBWFA_EXCITON_ANALYSIS_AD_H
#define LIBWFA_EXCITON_ANALYSIS_AD_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>
#include "exciton_moments.h"


namespace libwfa {

using namespace arma;

/** \brief Performs exciton analysis based on attachment and detachment density
        matrices.

    Performs an exciton analysis using the attachment and detachment density
    matrices.

    \ingroup libwfa
 **/
class exciton_analysis_ad {
private:
    const mom_builder_i &m_bld; //!< Builder for multipole moments
    const ab_matrix &m_adm; //!< Attachment density matrix
    const ab_matrix &m_ddm; //!< Detachment density matrix
    size_t m_mmax; //!< Max multipole moment to compute


public:
    exciton_analysis_ad(const mom_builder_i &bld,
        const ab_matrix &adm, const ab_matrix &ddm);

    /** \brief Perform analysis
        \param mom Resulting exciton moments
     **/
    void perform(ab_exciton_moments &mom) const;

private:
    void calculate(const arma::mat &adm,
        const arma::mat &ddm, exciton_moments &mom) const;
};

} // namespace libwfa


#endif
