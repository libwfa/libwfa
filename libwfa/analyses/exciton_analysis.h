#ifndef LIBWFA_EXCITON_ANALYSIS_H
#define LIBWFA_EXCITON_ANALYSIS_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>
#include "exciton_moments.h"


namespace libwfa {

/** \brief Computes the exciton moments based on a TDM

    Computes the exciton moments using the transition density matrix.

    \ingroup libwfa
 **/
class exciton_analysis {
private:
    exciton_moments *m_mom[2]; //!< Computed exciton moments

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
    ~exciton_analysis();

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
    static void calculate(const mom_builder_i &bld, const arma::mat &tdm,
        exciton_moments &mom);

    static void analysis(std::ostream &out, const exciton_moments &mom);

    static void print(std::ostream &out,
        const arma::vec &vec, size_t width = 10);
};

} // namespace libwfa


#endif
