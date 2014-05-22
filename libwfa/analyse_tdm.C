#include "analyse_tdm.h"

namespace libwfa {

using namespace arma;


analyse_tdm::analyse_tdm(const arma::Mat<double> &s, const ab_matrix &c,
    ev_printer_i &pr_nto) : m_s(s), m_nto(s, c, pr_nto) {

}


void analyse_tdm::do_register(const std::string &name,
    const ctnum_analysis_i &ana, const ctnum_printer_i &pr) {

    m_lst.insert(cna_map_t::value_type(name, cna(ana, pr)));
}


void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av,
    export_densities_i &dpr, export_orbitals_i &opr, std::ostream &out) {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, opr, std::cout);
    analyse_ctnum(tdm, out);

    dpr.perform(density_type::transition, tdm);
    dpr.perform(density_type::particle, eh.first);
    dpr.perform(density_type::hole, eh.second);

    av.first += eh.first;
    av.second += eh.second;
}


void analyse_tdm::analyse_ctnum(const ab_matrix &tdm, std::ostream &out) const {

    for (cna_map_t::const_iterator i = m_lst.begin(); i != m_lst.end(); i++) {

        out << i->first << std::endl;
        const cna &ctnum = i->second;

        ctnumbers(m_s, ctnum.analysis, ctnum.printer).perform(tdm, out);
    }
}


void analyse_tdm::analyse_nto(const ab_matrix &tdm, export_densities_i &dpr,
    export_orbitals_i &opr, std::ostream &out) const {

    ab_matrix_pair eh;
    m_nto.perform(tdm, eh, opr, out);

    dpr.perform(density_type::particle, eh.first);
    dpr.perform(density_type::hole, eh.second);
}


} // end namespace
