#ifndef LIBWFA_POP_ANALYSIS_AD_H
#define LIBWFA_POP_ANALYSIS_AD_H

#include <libwfa/core/ab_matrix.h>
#include "pop_analysis_i.h"
#include "pop_data.h"

namespace libwfa {

/** \brief Perform population analysis on state density matrix

    \ingroup libwfa
 **/
class pop_analysis_ad {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const pop_analysis_i &m_analysis;

public:
    /** \brief Constructor
        \param a Analysis object
     **/
    pop_analysis_ad(const pop_analysis_i &a) : m_analysis(a) { }

    /** \brief Perform population analysis
        \param[in] attach Attachment density matrix
        \param[in] detach Detachment density matrix
        \param[out] pop Population data
     **/
    void perform(const ab_matrix &attach, const ab_matrix &detach,
            pop_data &pr) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_SDM_H
