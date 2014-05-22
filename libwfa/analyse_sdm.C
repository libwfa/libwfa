#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include <libwfa/core/transformations_dm.h>
#include "analyse_sdm.h"

namespace libwfa {

using namespace arma;


analyse_sdm::analyse_sdm(const arma::Mat<double> &s,
    const ab_matrix &c, const ab_matrix &dm0, const pop_analysis_i &pop,
    export_densities_i &pr_d, export_orbitals_i &pr_o,
    ev_printer_i &pr_no, ev_printer_i &pr_ndo, pop_printer_i &pr_p) :
    m_dm0(dm0), m_ndo(s, c, pr_ndo), m_no(c, pr_no), m_pop(pop),
    m_pr_d(pr_d), m_pr_o(pr_o), m_pr_p(pr_p) {

}

void analyse_sdm::perform(const ab_matrix &dm, bool is_diff) {

    ab_matrix dm2(dm);
    if (is_diff) dm2 += m_dm0;
    else dm2 -= m_dm0;

    const ab_matrix &sdm(is_diff ? dm2 : dm);
    const ab_matrix &ddm(is_diff ? dm : dm2);

    analyse_no(sdm, m_pr_o, std::cout);

    ab_matrix_pair ad;
    analyse_ndo(ddm, ad, m_pr_o, std::cout);

    analyse_pop(sdm, ad, std::cout);

    m_pr_d.perform(density_type::state, sdm);
    m_pr_d.perform(density_type::attach, ad.first);
    m_pr_d.perform(density_type::detach, ad.second);
}


void analyse_sdm::analyse_pop(const ab_matrix &sdm,
        const ab_matrix_pair &ad, std::ostream &out) const {

    pop_data res;
    pop_analysis_dm(m_pop).perform(sdm, res);
    pop_analysis_ad(m_pop).perform(ad.first, ad.second, res);

    m_pr_p.perform(res, out);
}


} // end namespace
