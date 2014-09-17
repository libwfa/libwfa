#ifndef LIBWFA_EXCITON_ANALYSIS_H
#define LIBWFA_EXCITON_ANALYSIS_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>
#include "exciton_analysis_base.h"


namespace libwfa {

/** \brief Computes the exciton moments based on a TDM

    Computes the exciton moments using the transition density matrix.

    \ingroup libwfa
 **/
class exciton_analysis : public exciton_analysis_base {
public:
    /** \brief Constructor
        \param bld Moment builder
        \param tdm Transition density matrix
        \param maxmm Maximum moment (default: 2)
     **/
    exciton_analysis(const mom_builder_i &bld,
        const ab_matrix &tdm, size_t maxmm = 2);

    /** \brief Destructor
     **/
    virtual ~exciton_analysis() { }

    using exciton_analysis_base::analyse;

private:
    virtual void print_header(std::ostream &out, size_t off) const;

    virtual void analysis(std::ostream &out,
            const exciton_moments &mom, size_t off) const;

    static void calculate(const mom_builder_i &bld,
        const arma::mat &tdm, exciton_moments &mom);
};

} // namespace libwfa


#endif
