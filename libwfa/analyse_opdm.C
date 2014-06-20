#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include <libwfa/core/transformations_dm.h>
#include <memory>
#include "analyse_opdm.h"

namespace libwfa {

using namespace arma;

analyse_opdm::analyse_opdm(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &dm) : m_s(s), m_c(c), m_dm1(dm), m_dm2 (0),
            m_sdm (dm), m_ddm(dm){}

analyse_opdm::analyse_opdm(const arma::Mat<double> &s, const ab_matrix &c,
    const ab_matrix &dm0, const ab_matrix &dm, bool is_diff) :
    m_s(s), m_c(c), m_dm1(dm), m_dm2(build_dm(dm, dm0, is_diff)),
    m_sdm(is_diff ? *m_dm2 : m_dm1), m_ddm(is_diff ? m_dm1 : *m_dm2) {}


void analyse_opdm::do_register(const ev_printer_i &pr, bool is_no) {

    if (is_no) m_pr[0] = &pr;
    else if (m_dm2.get() != 0) m_pr[1] = &pr;
}


void analyse_opdm::do_register(const std::string &name,
    const pop_analysis_i &ana, const pop_printer_i &pr, pa_flag fl) {

    if (m_dm2.get() == 0) fl = (fl & pa_dm) == pa_dm ? pa_dm : pa_none;
    if (fl != pa_none) m_pa.insert(pa_map_t::value_type(name, pa(ana, pr, fl)));
}


void analyse_opdm::perform(export_data_i &pr, const contract_i &name,
        std::ostream &out) const {

    pr.perform(density_type::state, m_sdm);
    if (m_dm2.get() != 0) pr.perform(density_type::difference, m_ddm);

    if (m_pr[0] != 0) no_analysis(m_c, m_sdm, *m_pr[0]).perform(pr, out);

    ab_matrix at, de;
    if (m_pr[1] != 0) {
        ndo_analysis(m_s, m_c, m_ddm, *m_pr[1]).perform(at, de, pr, out);
        pr.perform(density_type::attach, at);
        pr.perform(density_type::detach, de);

        ex_analyse_ad analyse_ad;
        ex_ana_printer_ad ana_p_ad;

        analyse_ad.perform(at,de,name);
        ana_p_ad.perform(m_dm1.is_alpha_eq_beta(), analyse_ad, out);
    }

    for (pa_map_t::const_iterator i = m_pa.begin(); i != m_pa.end(); i++) {

        out << i->first << std::endl;
        const pa &pop = i->second;
        pop_data res;
        if ((pop.flag & pa_dm) == pa_dm)
            pop_analysis_dm(pop.analysis, m_sdm).perform(res);
        if (m_pr[1] != 0 && (pop.flag & pa_ad) == pa_ad)
            pop_analysis_ad(pop.analysis, at, de).perform(res);
        pop.printer.perform(res, out);
    }
}


std::auto_ptr<ab_matrix> analyse_opdm::build_dm(const ab_matrix &dm,
    const ab_matrix &dm0, bool is_diff) {

    std::auto_ptr<ab_matrix> dm2(new ab_matrix(dm));
    if (is_diff) *dm2 += dm0;
    else *dm2 -= dm0;
    return dm2;
}

} // end namespace
