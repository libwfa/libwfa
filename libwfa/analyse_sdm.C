#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include <libwfa/core/transformations_dm.h>
#include "analyse_sdm.h"

namespace libwfa {

using namespace arma;


analyse_sdm::analyse_sdm(const arma::Mat<double> &s,
    const ab_matrix &c, const ab_matrix &dm0,
    const ev_printer_i &pr_no, const ev_printer_i &pr_ndo) :
    m_dm0(dm0), m_no(c, pr_no), m_ndo(s, c, pr_ndo) {

}


void analyse_sdm::do_register(const std::string &name,
    const pop_analysis_i &ana, const pop_printer_i &pr, pa_flag fl) {

    m_lst.insert(pa_map_t::value_type(name, pa(ana, pr, fl)));
}


void analyse_sdm::perform(const ab_matrix &dm, export_densities_i &dpr,
    export_orbitals_i &opr, std::ostream &out, bool is_diff) {

    ab_matrix dm2(dm);
    if (is_diff) dm2 += m_dm0;
    else dm2 -= m_dm0;

    const ab_matrix &sdm(is_diff ? dm2 : dm);
    const ab_matrix &ddm(is_diff ? dm : dm2);

    analyse_no(sdm, opr, out);

    ab_matrix_pair ad;
    analyse_ndo(ddm, ad, opr, out);

    analyse_pop(sdm, ad, out);

    dpr.perform(density_type::state, sdm);
    dpr.perform(density_type::attach, ad.first);
    dpr.perform(density_type::detach, ad.second);
}


void analyse_sdm::analyse_pop(const ab_matrix &sdm,
    const ab_matrix_pair &ad, std::ostream &out) const {

    for (pa_map_t::const_iterator i = m_lst.begin(); i != m_lst.end(); i++) {

        out << i->first << std::endl;
        const pa &pop = i->second;
        pop_data res;
        if ((pop.flag & pa_dm) == pa_dm) {
            pop_analysis_dm(pop.analysis).perform(sdm, res);
        }
        if ((i->second.flag & pa_ad) == pa_ad) {
            pop_analysis_ad(pop.analysis).perform(ad.first, ad.second, res);
        }
        pop.printer.perform(res, out);
    }
}


} // end namespace
