#include "analyse_tdm.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av,
    export_densities_i &dm_print, export_orbitals_i &nto_print,
    ev_data_i &prn, ctnum_data_i &prct) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, nto_print, prn);
    m_ct.perform(tdm, prct);

    dm_print.perform(density_type::transition, tdm);
    dm_print.perform(density_type::particle, eh.first);
    dm_print.perform(density_type::hole, eh.second);

    av.first += eh.first;
    av.second += eh.second;
}


void analyse_tdm::nto_analysis(const ab_matrix &tdm,
    export_densities_i &dm_print, export_orbitals_i &nto_print,
    ev_data_i &pr) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, nto_print, pr);

    dm_print.perform(density_type::particle, eh.first);
    dm_print.perform(density_type::hole, eh.second);
}


} // end namespace
