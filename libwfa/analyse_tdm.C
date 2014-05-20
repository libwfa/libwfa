#include "analyse_tdm.h"

namespace libwfa {

using namespace arma;


void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av,
    export_densities_i &pr_d, export_orbitals_i &pr_o,
    ev_printer_i &pr_e, ctnum_printer_i &pr_c) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, pr_o, pr_e);
    m_ct.perform(tdm, pr_c);

    pr_d.perform(density_type::transition, tdm);
    pr_d.perform(density_type::particle, eh.first);
    pr_d.perform(density_type::hole, eh.second);

    av.first += eh.first;
    av.second += eh.second;
}


void analyse_tdm::nto_analysis(const ab_matrix &tdm,
    export_densities_i &pr_d, export_orbitals_i &pr_o,
    ev_printer_i &pr_e) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, pr_o, pr_e);

    pr_d.perform(density_type::particle, eh.first);
    pr_d.perform(density_type::hole, eh.second);
}


} // end namespace
