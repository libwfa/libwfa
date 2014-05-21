#include "analyse_tdm.h"

namespace libwfa {

using namespace arma;


void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av) {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, m_pr_o, m_pr_nto);
    m_ct.perform(tdm, m_pr_ct);

    m_pr_d.perform(density_type::transition, tdm);
    m_pr_d.perform(density_type::particle, eh.first);
    m_pr_d.perform(density_type::hole, eh.second);

    av.first += eh.first;
    av.second += eh.second;
}


void analyse_tdm::analyse_nto(const ab_matrix &tdm,
    export_densities_i &pr_d, export_orbitals_i &pr_o,
    ev_printer_i &pr_nto) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, pr_o, pr_nto);

    pr_d.perform(density_type::particle, eh.first);
    pr_d.perform(density_type::hole, eh.second);
}


} // end namespace
