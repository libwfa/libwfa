#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include "analyse_ddm.h"

namespace libwfa {

using namespace arma;


analyse_ddm::analyse_ddm(const arma::Mat<double> &s, const ab_matrix &c,
    const ab_matrix &ddm, const ev_printer_i &prndo) :
    m_s(s), m_c(c), m_ddm(ddm), m_prndo(prndo) {

}


void analyse_ddm::do_register(const std::string &name,
    const pop_analysis_i &ana, const pop_printer_i &pr) {

    m_lst.insert(pa_map_t::value_type(name, pa(ana, pr)));
}


void analyse_ddm::perform(export_data_i &pr, std::ostream &out) const {

    ab_matrix at, de;
    ndo_analysis(m_s, m_c, m_ddm, m_prndo).perform(at, de, pr, out);

    for (pa_map_t::const_iterator i = m_lst.begin(); i != m_lst.end(); i++) {

        out << i->first << std::endl;
        const pa &pop = i->second;
        pop_data res;
        pop_analysis_ad(pop.analysis, at, de).perform(res);
        pop.printer.perform(res, out);
    }

    pr.perform(density_type::difference, m_ddm);
    pr.perform(density_type::attach, at);
    pr.perform(density_type::detach, de);
}


} // end namespace
