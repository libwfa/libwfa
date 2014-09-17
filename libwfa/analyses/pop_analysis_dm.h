#ifndef LIBWFA_POP_ANALYSIS_DM_H
#define LIBWFA_POP_ANALYSIS_DM_H

#include <libwfa/core/ab_matrix.h>
#include "pop_analysis_i.h"
#include "pop_data.h"

namespace libwfa {

/** \brief Perform population analysis on state density matrix

    Performs a population analysis of the state density matrix provided.
    The analysis is performed on the spin averaged density matrix and,
    if unrestricted, on the spin difference density matrix.
    The base population (e.g. nuclear charges) is added to the spin average
    density matrix.

    \ingroup libwfa
 **/
class pop_analysis_dm {
private:
    const pop_analysis_i &m_analysis; //!< Analysis
    const arma::vec &m_p0; //!< Base population
    const ab_matrix &m_sdm; //!< State density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param sdm State density matrix
        \param p0 Base population
     **/
    pop_analysis_dm(const pop_analysis_i &a, const arma::vec &p0,
        const ab_matrix &sdm) : m_analysis(a), m_p0(p0), m_sdm(sdm) { }

    /** \brief Perform population analysis
        \param pop Resulting population data
     **/
    void perform(pop_data &pop) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_DM_H
