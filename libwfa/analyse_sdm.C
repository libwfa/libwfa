#include "analyse_sdm.h"
#include "pop_analysis_ad.h"
#include "pop_analysis_dm.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void analyse_sdm::perform(const ab_matrix &sdm, export_densities_i &pr_d,
    export_orbitals_i &pr_o, ev_printer_i &pr_e1, ev_printer_i &pr_e2,
    pop_printer_i &pr_p) const {

    no_analysis(sdm, pr_o, pr_e1);

    ab_matrix_pair ad;
    ndo_analysis(sdm, ad, pr_o, pr_e2);

    pop_analysis(sdm, ad, pr_p);

    pr_d.perform(density_type::state, sdm);
    pr_d.perform(density_type::attach, ad.first);
    pr_d.perform(density_type::detach, ad.second);
}


void analyse_sdm::pop_analysis(const ab_matrix &sdm,
        const ab_matrix_pair &ad, pop_printer_i &pr) const {

    pop_data res;
    pop_analysis_dm(m_pop).perform(sdm, res);
    pop_analysis_ad(m_pop).perform(ad.first, ad.second, res);

    pr.perform(res);
}


} // end namespace
