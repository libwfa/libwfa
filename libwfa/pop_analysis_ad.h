#ifndef LIBWFA_POP_ANALYSIS_AD_H
#define LIBWFA_POP_ANALYSIS_AD_H

#include "ab_matrix.h"
#include "pop_analysis_i.h"
#include "pop_print_i.h"

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
        \param[in] ad Attachment and detachment density matrix (attach first)
        \param[out] pr Printer of population data
     **/
    void perform(const ab_matrix_pair &ad, pop_print_i &pr) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_SDM_H
