#include "analyse_tdm.h"

namespace libwfa {

using namespace arma;


analyse_tdm::analyse_tdm(const arma::Mat<double> &s,
    const ab_matrix &c, const ctnum_analysis_i &ctnum,
    export_densities_i &pr_d, export_orbitals_i &pr_o,
    ev_printer_i &pr_nto, ctnum_printer_i &pr_ct) :
    m_ct(s, ctnum, pr_ct), m_nto(s, c, pr_nto), m_pr_d(pr_d), m_pr_o(pr_o) {

}

void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av) {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, m_pr_o, std::cout);
    m_ct.perform(tdm, std::cout);

    m_pr_d.perform(density_type::transition, tdm);
    m_pr_d.perform(density_type::particle, eh.first);
    m_pr_d.perform(density_type::hole, eh.second);

    av.first += eh.first;
    av.second += eh.second;
}


void analyse_tdm::analyse_nto(const ab_matrix &tdm, export_densities_i &pr_d,
    export_orbitals_i &pr_o, std::ostream &out) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, pr_o, out);

    pr_d.perform(density_type::particle, eh.first);
    pr_d.perform(density_type::hole, eh.second);
}


} // end namespace
