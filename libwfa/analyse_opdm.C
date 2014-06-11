#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include <libwfa/core/transformations_dm.h>
#include "analyse_opdm.h"

namespace libwfa {

using namespace arma;


analyse_opdm::analyse_opdm(const arma::Mat<double> &s, const ab_matrix &c,
    const ab_matrix &dm0, const ab_matrix &dm, bool is_diff) :
    m_s(s), m_c(c), m_dm0(dm0), m_dm(dm), m_no(0), m_ndo(0), m_is_diff(is_diff) {

}


void analyse_opdm::do_register(orbital_type ot, const ev_printer_i &pr) {

    if (ot == orbital_type::no) m_no = &pr;
    else if (ot == orbital_type::ndo) m_ndo = &pr;
}


void analyse_opdm::do_register(const std::string &name,
    const pop_analysis_i &ana, const pop_printer_i &pr, pa_flag fl) {

    m_pa.insert(pa_map_t::value_type(name, pa(ana, pr, fl)));
}


void analyse_opdm::perform(export_data_i &pr, std::ostream &out) const {

    ab_matrix dm2(m_dm);
    if (m_is_diff) dm2 += m_dm0;
    else dm2 -= m_dm0;

    const ab_matrix &sdm(m_is_diff ? dm2 : m_dm);
    const ab_matrix &ddm(m_is_diff ? m_dm : dm2);

    pr.perform(density_type::state, sdm);
    pr.perform(density_type::difference, ddm);

    if (m_no != 0) {
        no_analysis(m_c, sdm, *m_no).perform(pr, out);
    }

    ab_matrix at, de;
    if (m_ndo != 0) {
        ndo_analysis(m_s, m_c, ddm, *m_ndo).perform(at, de, pr, out);
        pr.perform(density_type::attach, at);
        pr.perform(density_type::detach, de);
    }

    for (pa_map_t::const_iterator i = m_pa.begin(); i != m_pa.end(); i++) {

        out << i->first << std::endl;
        const pa &pop = i->second;
        pop_data res;
        if ((pop.flag & pa_dm) == pa_dm)
            pop_analysis_dm(pop.analysis, sdm).perform(res);
        if (m_ndo != 0 && (pop.flag & pa_ad) == pa_ad)
            pop_analysis_ad(pop.analysis, at, de).perform(res);
        pop.printer.perform(res, out);
    }
}


} // end namespace
