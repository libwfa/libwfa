#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

#include "ctnumbers.h"
#include "nto_analysis.h"

namespace libwfa {

/** \brief Class to perform transition density matrix analysis

    Performs the analyses of the current transition density according to the
    parameters.

    \ingroup libwfa
 **/
class analyse_tdm {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    ctnumbers m_ct;
    nto_analysis m_nto;

public:
    analyse_tdm(const ctnum_analysis_i &ctnum,
        const arma::Mat<double> &s, const ab_matrix &c) :
        m_ct(ctnum, s), m_nto(s, c) { }


    void perform(const ab_matrix &tdm, ab_matrix_pair &av, std::ostream &out,
        export_densities_i &dm_print, export_orbitals_i &nto_print,
        nto_data_i &prn, ctnum_data_i &prct);
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
