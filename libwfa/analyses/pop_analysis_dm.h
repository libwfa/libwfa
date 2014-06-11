#ifndef LIBWFA_POP_ANALYSIS_DM_H
#define LIBWFA_POP_ANALYSIS_DM_H

#include <libwfa/core/ab_matrix.h>
#include "pop_analysis_i.h"
#include "pop_data.h"

namespace libwfa {

/** \brief Perform population analysis on state density matrix

    \ingroup libwfa
 **/
class pop_analysis_dm {
private:
    const pop_analysis_i &m_analysis; //!< Analysis
    const ab_matrix &m_sdm; //!< State density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param sdm State density matrix
     **/
    pop_analysis_dm(const pop_analysis_i &a, const ab_matrix &sdm) :
        m_analysis(a), m_sdm(sdm) { }

    /** \brief Perform population analysis
        \param pop Resulting population data
     **/
    void perform(pop_data &pop) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_DM_H
