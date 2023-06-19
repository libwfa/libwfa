#ifndef LIBWFA_POP_ANALYSIS_TDM_H
#define LIBWFA_POP_ANALYSIS_TDM_H

#include <libwfa/core/ab_matrix.h>
#include "pop_analysis_i.h"
#include "pop_data.h"

namespace libwfa {

/** \brief Perform population analysis on transition density matrix

    Performs a population analysis of the transition density matrix provided.
    This computes transition charges of the spin-averaged DM.

    \ingroup libwfa
 **/
class pop_analysis_tdm {
private:
    const pop_analysis_i &m_analysis; //!< Analysis
    const ab_matrix &m_tdm; //!< Transition density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param tdm State density matrix
     **/
    pop_analysis_tdm(const pop_analysis_i &a,
        const ab_matrix &tdm) : m_analysis(a), m_tdm(tdm) { }

    /** \brief Perform population analysis
        \param pop Resulting population data
     **/
    void perform(pop_data &pop) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_TDM_H
