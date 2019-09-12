#ifndef LIBWFA_EXCITON_ANALYSIS_AD_H
#define LIBWFA_EXCITON_ANALYSIS_AD_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>
#include "exciton_analysis_base.h"


namespace libwfa {

using namespace arma;

/** \brief Performs exciton analysis based on attachment and detachment density
        matrices.

    Performs an exciton analysis using the attachment and detachment density
    matrices.

    \ingroup libwfa
 **/
class exciton_analysis_ad : public exciton_analysis_base {
public:
    /** \brief Constructor
        \param bld Moment builder
        \param adm Attachment density matrix
        \param ddm Detachment density matrix
        \param maxmm Highest moment to perform analysis
     **/
    exciton_analysis_ad(const mom_builder_i &bld,
        const ab_matrix &adm, const ab_matrix &ddm, size_t maxmm = 2);

    /** \brief Destructor
     **/
    virtual ~exciton_analysis_ad() { }

    using exciton_analysis_base::analyse;

    virtual void combine(const exciton_moments &a, const exciton_moments &b,
        exciton_moments &res) const;

private:
    virtual void print_header(std::ostream &out, size_t off) const;

    virtual void analysis(std::ostream &out,
            const exciton_moments &mom, size_t off) const;

    static void calculate(const mom_builder_i &bld, const arma::mat &adm,
        const arma::mat &ddm, exciton_moments &mom);
};

} // namespace libwfa


#endif
