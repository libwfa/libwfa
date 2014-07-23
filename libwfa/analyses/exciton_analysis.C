#include "exciton_analysis.h"

namespace libwfa {

using namespace arma;


exciton_analysis::exciton_analysis(const mom_builder_i &bld,
    const ab_matrix &tdm) : m_bld(bld), m_tdm(tdm), m_mmax(2) {

    m_mmax = std::min(m_bld.max_moment(), m_mmax);
}


void exciton_analysis::perform(ab_exciton_moments &mom) const {

    bool aeqb = m_tdm.is_alpha_eq_beta();
    if (aeqb) mom.set_alpha_eq_beta();
    else mom.set_alpha_neq_beta();

    calculate(m_tdm.alpha(), mom.alpha());
    if (! aeqb) {
        calculate(m_tdm.beta(), mom.beta());
    }
}


void exciton_analysis::calculate(const arma::Mat<double> &tdm,
    exciton_moments &mom) const {

    double n = m_bld.perform(tdm, 0, 0, 0, 0);

    mom = exciton_moments(m_mmax);
    for (size_t i = 1; i <= m_mmax; i++) {

        for (size_t j = 0; j <= i; j++) {

            Col<double> mj(3, fill::zeros);
            for (size_t k = 0; k < 3; k++)
                mj(k) = m_bld.perform(tdm, k, j, k, i - j) / n;
            mom.set(j, i - j, mj);
        }
    }
}


}// end namespace libwfa





