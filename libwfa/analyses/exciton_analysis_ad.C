#include "exciton_analysis_ad.h"

namespace libwfa {

using namespace arma;


exciton_analysis_ad::exciton_analysis_ad(const mom_builder_i &bld,
    const ab_matrix &adm, const ab_matrix &ddm) :
    m_bld(bld), m_adm(adm), m_ddm(ddm), m_mmax(2) {

    m_mmax = std::min(m_bld.max_moment(), m_mmax);
}

void exciton_analysis_ad::perform(ab_exciton_moments &mom) const {

    bool aeqb = m_adm.is_alpha_eq_beta();
    if (aeqb) mom.set_alpha_eq_beta();
    else mom.set_alpha_neq_beta();

    calculate(m_adm.alpha(), m_ddm.alpha(), mom.alpha());
    if (! aeqb) calculate(m_adm.beta(), m_ddm.beta(), mom.beta());
}


void exciton_analysis_ad::calculate(const arma::mat &adm,
    const arma::mat &ddm, exciton_moments &mom) const {

    double n = m_bld.perform(adm, 'x', 0);

    mom = exciton_moments(m_mmax);
    for (size_t i = 1; i <= m_mmax; i++) {
        vec mh(3, fill::zeros), me(3, fill::zeros);
        for (size_t k = 0; k < 3; k++) {
            me(k) = m_bld.perform( adm, k, i) / n;
            mh(k) = m_bld.perform(-ddm, k, i) / n;
        }
        mom.set(i, 0, me);
        mom.set(0, i, mh);
    }
}


}// end namespace libwfa





