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
    exciton_moments *m_mom[2]; //!< Computed exciton moments

public:
    exciton_analysis_ad(const mom_builder_i &bld,
        const ab_matrix &adm, const ab_matrix &ddm, size_t maxmm = 2);

    /** \brief Destructor
     **/
    ~exciton_analysis_ad();

    /** \brief Computed exciton moment
        \param spin If true: beta spin; else alpha spin
     **/
    const exciton_moments &moment(bool spin) {
        return *m_mom[(spin && m_mom[1] ? 1 : 0)];
    }

    /** \brief Perform analysis
        \param out Output stream
     **/
    void analyse(std::ostream &out) const;

private:
    static void calculate(const mom_builder_i &bld, const arma::mat &adm,
        const arma::mat &ddm, exciton_moments &mom);

    static void analysis(std::ostream &out, const exciton_moments &mom);

    static void print(std::ostream &out,
        const arma::vec &vec, size_t width = 10);
};

} // namespace libwfa


#endif
