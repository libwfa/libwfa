#ifndef LIBWFA_POP_ANALYSIS_AD_H
#define LIBWFA_POP_ANALYSIS_AD_H

#include <libwfa/core/ab_matrix.h>
#include "pop_analysis_i.h"
#include "pop_data.h"

namespace libwfa {

/** \brief Perform population analysis on a/d density matrices

    \ingroup libwfa
 **/
class pop_analysis_ad {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const pop_analysis_i &m_analysis; //!< Analysis
    const ab_matrix &m_at; //!< Attachment density matrix
    const ab_matrix &m_de; //!< Detachment density matrix

public:
    /** \brief Constructor
        \param a Analysis object
        \param at Attachment density matrix
        \param de Detachment density matrix
     **/
    pop_analysis_ad(const pop_analysis_i &a, const ab_matrix &at,
        const ab_matrix &de) : m_analysis(a), m_at(at), m_de(de) { }

    /** \brief Perform population analysis
        \param pop Population data
     **/
    void perform(pop_data &pr) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_SDM_H
