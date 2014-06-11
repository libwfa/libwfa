#include <libwfa/analyses/ctnumbers.h>
#include <libwfa/analyses/nto_analysis.h>
#include "analyse_tdm.h"

namespace libwfa {

using namespace arma;


analyse_tdm::analyse_tdm(const Mat<double> &s, const ab_matrix &c,
    const ab_matrix &tdm, const ev_printer_i &prnto) :
    m_s(s), m_c(c), m_tdm(tdm), m_prnto(prnto) {

}


void analyse_tdm::do_register(const std::string &name,
    const ctnum_analysis_i &ana, const ctnum_printer_i &pr) {

    m_lst.insert(cna_map_t::value_type(name, cna(ana, pr)));
}


void analyse_tdm::perform(ab_matrix &edm_av, ab_matrix &hdm_av,
    export_data_i &pr, std::ostream &out) {

    {
    ab_matrix edm, hdm;
    nto_analysis(m_s, m_c, m_tdm, m_prnto).perform(edm, hdm, pr, out);

    pr.perform(density_type::transition, m_tdm);
    pr.perform(density_type::particle, edm);
    pr.perform(density_type::hole, hdm);

    edm_av += edm;
    hdm_av += hdm;
    }

    for (cna_map_t::const_iterator i = m_lst.begin(); i != m_lst.end(); i++) {

        out << i->first << std::endl;
        const cna &ctnum = i->second;

        ab_matrix om;
        double om_tot[2];
        ctnumbers(ctnum.analysis, m_s, m_tdm).perform(om, om_tot);
        ctnum.printer.perform(om, om_tot, out);
    }
}


} // end namespace
