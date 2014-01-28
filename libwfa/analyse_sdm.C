#include "analyse_sdm.h"
#include "pop_analysis_ad.h"
#include "pop_analysis_dm.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void analyse_sdm::perform(const ab_matrix &sdm, export_densities_i &dm_print,
    export_orbitals_i &orb_print, ev_data_i &no, ev_data_i &ndo,
    pop_print_i &pr) const {

    no_analysis(sdm, orb_print, no);

    ab_matrix_pair ad;
    ndo_analysis(sdm, ad, orb_print, ndo);

    pop_analysis(sdm, ad, pr);

    dm_print.perform(density_type::state, sdm);
    dm_print.perform(density_type::attach, ad.first);
    dm_print.perform(density_type::detach, ad.second);
}


void analyse_sdm::pop_analysis(const ab_matrix &sdm,
        const ab_matrix_pair &ad, pop_print_i &pr) const {

    pop_data res;
    pop_analysis_dm(m_pop).perform(sdm, res);
    pop_analysis_ad(m_pop).perform(ad.first, ad.second, res);

    pr.perform(res);
}


} // end namespace
