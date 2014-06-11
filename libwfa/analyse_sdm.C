#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include <libwfa/core/transformations_dm.h>
#include "analyse_sdm.h"

namespace libwfa {

using namespace arma;


analyse_sdm::analyse_sdm(const ab_matrix &c,
    const ab_matrix &sdm, const ev_printer_i &prno) :
    m_c(c), m_sdm(sdm), m_prno(prno) {

}


void analyse_sdm::do_register(const std::string &name,
    const pop_analysis_i &ana, const pop_printer_i &pr) {

    m_lst.insert(pa_map_t::value_type(name, pa(ana, pr)));
}


void analyse_sdm::perform(export_data_i &pr, std::ostream &out) const {

    no_analysis(m_c, m_sdm, m_prno).perform(pr, out);

    for (pa_map_t::const_iterator i = m_lst.begin(); i != m_lst.end(); i++) {

        out << i->first << std::endl;
        const pa &pop = i->second;
        pop_data res;
        if ((pop.flag & pa_dm) == pa_dm) {
            pop_analysis_dm(pop.analysis).perform(sdm, res);
        }
        if ((pop.flag & pa_ad) == pa_ad) {
            pop_analysis_ad(pop.analysis).perform(ad.first, ad.second, res);
        }
        pop.printer.perform(res, out);
    }

    pr.perform(density_type::state, m_sdm);
}


} // end namespace
