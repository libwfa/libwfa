#include <libwfa/analyses/ctnumbers.h>
#include <libwfa/analyses/exciton_analysis_ad.h>
#include <libwfa/analyses/exciton_analysis.h>
#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/nto_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include "wf_analysis.h"

namespace libwfa {

using namespace arma;


wf_analysis::wf_analysis(wf_analysis_data_i *h,
    const std::set<std::string> &p1, const wfa_params &p2) :
    m_h(h), m_p2(p2), m_init_av(false), m_sa(0) {

    for (std::set<std::string>::const_iterator i = p1.begin();
            i != p1.end(); i++) {
        const std::string &a = *i; activate(a);
    }
}


void wf_analysis::activate(const std::string &a) {

    ana_type p = convert(a);
    if (p == NA) return;

    m_p1.set(p);
    if (p == FORM_AD)
        m_p1.set(NDO);
    else if (p == EXCITON_AD) {
        m_p1.set(NDO);
        m_p1.set(FORM_AD);
    }
    else if (p == SA_NTO)
        m_p1.set(FORM_EH);
}


void wf_analysis::deactivate(const std::string &a) {

    ana_type p = convert(a);
    if (p == NA) return;

    m_p1.reset(p);
    if (p == NDO) {
        m_p1.reset(EXCITON_AD);
        m_p1.reset(FORM_AD);
    }
    else if (p == FORM_AD)
        m_p1.reset(EXCITON_AD);
    else if (p == FORM_EH)
        m_p1.reset(SA_NTO);
}


void wf_analysis::analyse_opdm(std::ostream &out, const std::string &name,
    const std::string &desc, const ab_matrix &ddm, const ab_matrix &dm0) {

    ab_matrix sdm(ddm);
    sdm += dm0;

    // Create printer for orbitals and densities
    std::auto_ptr<density_printer_i> pr1(m_h->density_printer(name, desc));
    std::auto_ptr<orbital_printer_i> pr2(m_h->orbital_printer(name, desc));

    // Export density matrices first
    pr1->perform(density_type::state, sdm);
    pr1->perform(density_type::difference, ddm);

    // Get prerequisites of further analyses
    const arma::mat &s = m_h->overlap();
    const ab_matrix &c = m_h->coefficients();

    if (m_p1.test(NO)) {
        no_analysis no(s, c, sdm);
        no.analyse(out, m_p2.nno);
        no.export_orbitals(*pr2);
        out << std::endl;
    }

    ab_matrix at, de;
    if (m_p1.test(NDO)) {
        ndo_analysis ndo(s, c, ddm);
        ndo.analyse(out, m_p2.nndo);
        ndo.export_orbitals(*pr2, m_p2.nndo);
        out << std::endl;

        if (m_p1.test(FORM_AD)) {
            ndo.form_ad(at, de);
            pr1->perform(density_type::attach, at);
            pr1->perform(density_type::detach, de);
        }
    }

    // Perform population analyses
    for (size_t i = 0; i < m_h->n_pop_analyses(); i++) {

        pop_data pdata;
        const pop_analysis_i &pa = m_h->pop_analysis(i);
        const vec &p0 = m_h->ref_population(i);
        const std::string &pname = m_h->pop_name(i);
        const std::vector<std::string> &l = m_h->pop_labels(i);

        pop_analysis_dm(pa, p0, sdm).perform(pdata);
        if (m_p1.test(FORM_AD))
            pop_analysis_ad(pa, at, de).perform(pdata);

        out << pname << std::endl;
        pdata.print(out, l);
        out << std::endl;
    }

    if (m_p1.test(EXCITON_AD)) {
        exciton_analysis_ad(m_h->mom_builder(), at, de).analyse(out, 0);
        out << std::endl;
    }
}


void wf_analysis::analyse_opdm(std::ostream &out, const std::string &name,
    const std::string &desc, const ab_matrix &sdm) {

    wf_analysis_data_i &h = *wf_analysis::m_h;

    // Create printer for orbitals and densities
    std::auto_ptr<density_printer_i> pr1(h.density_printer(name, desc));
    std::auto_ptr<orbital_printer_i> pr2(h.orbital_printer(name, desc));

    // Export density matrices first
    pr1->perform(density_type::state, sdm);

    // Get prerequisites of further analyses
    const arma::mat &s = h.overlap();
    const ab_matrix &c = h.coefficients();

    //opdm_params p1 = m_h->get_opdm_params();
    if (m_p1.test(NO)) {
        no_analysis no(s, c, sdm);
        no.analyse(out, m_p2.nno);
        no.export_orbitals(*pr2);
        out << std::endl;
    }

    // Perform population analyses
    for (size_t i = 0; i < h.n_pop_analyses(); i++) {

        pop_data pdata;
        const pop_analysis_i &pa = h.pop_analysis(i);
        const vec &p0 = h.ref_population(i);
        const std::string &pname = h.pop_name(i);
        const std::vector<std::string> &l = h.pop_labels(i);

        pop_analysis_dm(pa, p0, sdm).perform(pdata);
        out << pname << std::endl;
        pdata.print(out, l);
        out << std::endl;
    }

}


void wf_analysis::analyse_optdm(std::ostream &out, const std::string &name,
    const std::string &desc, const ab_matrix &tdm) {

    // Create printer for orbitals and densities
    std::auto_ptr<density_printer_i> pr1(m_h->density_printer(name, desc));
    std::auto_ptr<orbital_printer_i> pr2(m_h->orbital_printer(name, desc));

    // Export density matrices first
    pr1->perform(density_type::transition, tdm);

    // Get prerequisites of further analyses
    const arma::mat &s = m_h->overlap();
    const ab_matrix &c = m_h->coefficients();

    // If NTO formatter exists, do NTO analysis
    if (m_p1.test(FORM_EH)) {

        ab_matrix edm, hdm;
        nto_analysis::form_eh(s, tdm, edm, hdm);

        if (m_p1.test(NTO)) {
            nto_analysis nto(s, c, edm, hdm);
            nto.analyse(out, m_p2.nnto);
            nto.export_orbitals(*pr2, m_p2.nto_thresh);
            out << std::endl;
        }
        pr1->perform(density_type::particle, edm);
        pr1->perform(density_type::hole, hdm);
        if (m_p1.test(SA_NTO)) add_to_average(edm, hdm);
    }
    else if (m_p1.test(NTO)) {

        nto_analysis nto(s, c, tdm);
        nto.analyse(out, m_p2.nnto);
        nto.export_orbitals(*pr2, m_p2.nto_thresh);
        out << std::endl;
    }


    if (m_h->n_ctnum_analyses() != 0) {

        for (size_t i = 0; i < m_h->n_ctnum_analyses(); i++) {

            const ctnum_analysis_i &ca = m_h->ctnum_analysis(i);
            const std::string &cname = m_h->ctnum_name(i);
            std::auto_ptr<ctnum_printer_i> cpr(m_h->ctnum_printer(i, name, desc));

            out << cname << std::endl;
            ctnumbers ct(ca, tdm);
            ct.analyse(out);
            ct.do_export(*cpr);
            out << std::endl;
        }
    }

    if (m_p1.test(EXCITON)) {
        exciton_analysis(m_h->mom_builder(), tdm).analyse(out, 0);
        out << std::endl;
    }
}


bool wf_analysis::setup_sa_ntos(std::ostream &out) {

    if (! m_p1.test(SA_NTO) || ! m_init_av) return false;
    if (m_sa.get()) return true;

    ab_matrix edm(true), hdm(true);
    edm.alpha() = 0.5 * (m_edm_av.alpha() + m_edm_av.beta());
    hdm.alpha() = 0.5 * (m_hdm_av.alpha() + m_hdm_av.beta());

    const arma::mat &s = m_h->overlap();
    const ab_matrix &c = m_h->coefficients();

    nto_analysis sa_ntos(s, c, edm, hdm);
    if (m_p1.test(NTO)) {
        std::auto_ptr<orbital_printer_i> pr(m_h->orbital_printer("sa_nto",
                "State-averaged NTOs"));
        sa_ntos.analyse(out, m_p2.nnto);
        sa_ntos.export_orbitals(*pr, m_p2.nto_thresh);
    }
    m_sa = std::auto_ptr<sa_nto_analysis>(new sa_nto_analysis(s, sa_ntos));

    return true;
}

bool wf_analysis::post_process_optdm(std::ostream &out, const ab_matrix &tdm) {

    if (m_sa.get() == 0) return false;

    out << "Decomposition into state-averaged NTOs" << std::endl;
    m_sa->analyse(out, tdm);
    return true;
}


void wf_analysis::add_to_average(const ab_matrix &edm, const ab_matrix &hdm) {

    if (m_init_av) {
        m_edm_av += edm;
        m_hdm_av += hdm;
    }
    else {
        m_edm_av = edm;
        m_hdm_av = hdm;
        m_init_av = true;
    }
}

wf_analysis::ana_type wf_analysis::convert(const std::string &a) {

    // TODO: Improve conversion
    if (a == "no") return NO;
    else if (a == "ndo") return NDO;
    else if (a == "ad") return FORM_AD;
    else if (a == "exciton_ad") return EXCITON_AD;
    else if (a == "nto") return NTO;
    else if (a == "eh") return FORM_EH;
    else if (a == "exciton") return EXCITON;
    else if (a == "sa_nto") return SA_NTO;
    else return NA;
}

std::auto_ptr<wf_analysis> wf_analysis_static::analysis(0);

} // namespace libwfa
