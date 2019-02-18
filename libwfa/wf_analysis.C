#include <libwfa/analyses/ctnumbers.h>
#include <libwfa/analyses/dens_mom.h>
#include <libwfa/analyses/exciton_analysis_ad.h>
#include <libwfa/analyses/exciton_analysis.h>
#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/nto_analysis.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_analysis_dm.h>
#include "wf_analysis.h"
#include <fstream>
#include <iomanip>

namespace libwfa {

using namespace arma;


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

    const arma::mat &s = m_h->overlap();
    const ab_matrix &c = m_h->coefficients();

    // Get prerequisites of further analyses
    if (m_h->is_active(wf_analysis_data_i::NO)) {
        orbital_params pno = m_h->get_orbital_params(orbital_type::NO);
        no_analysis no(s, c, sdm);
        no.analyse(out, pno.norb);
        no.export_orbitals(*pr2, pno.thresh);
        out << std::endl;
    }

    ab_matrix at, de;
    if (m_h->is_active(wf_analysis_data_i::NDO)) {
        orbital_params pndo = m_h->get_orbital_params(orbital_type::NDO);
        ndo_analysis ndo(s, c, ddm);
        ndo.analyse(out, pndo.norb);
        ndo.export_orbitals(*pr2, pndo.thresh);
        out << std::endl;

        if (m_h->is_active(wf_analysis_data_i::FORM_AD)) {
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
        if (m_h->is_active(wf_analysis_data_i::FORM_AD))
            pop_analysis_ad(pa, at, de).perform(pdata);

        out << pname << std::endl;
        pdata.print(out, l);
        out << std::endl;
    }

    if (m_h->is_active(wf_analysis_data_i::DENS_MOM)) {
        dens_mom(m_h->mom_builder(), sdm, m_h->coordinates(), m_h->atomic_charges()).analyse(out, 0);
    }

    if (m_h->is_active(wf_analysis_data_i::EXCITON_AD)) {
        exciton_analysis_ad(m_h->mom_builder(), at, de).analyse(out, 0);
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
    if (m_h->is_active(wf_analysis_data_i::NO)) {
        orbital_params pno = m_h->get_orbital_params(orbital_type::NO);
        no_analysis no(s, c, sdm);
        no.analyse(out, pno.norb);
        no.export_orbitals(*pr2, pno.thresh);
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

    if (m_h->is_active(wf_analysis_data_i::DENS_MOM)) {
        dens_mom(m_h->mom_builder(), sdm, m_h->coordinates(), m_h->atomic_charges()).analyse(out, 0);
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
    if (m_h->is_active(wf_analysis_data_i::FORM_EH)) {

        ab_matrix edm, hdm;
        nto_analysis::form_eh(s, tdm, edm, hdm);

        if (m_h->is_active(wf_analysis_data_i::NTO)) {
            orbital_params pnto = m_h->get_orbital_params(orbital_type::NTO);
            nto_analysis nto(s, c, edm, hdm);
            nto.analyse(out, pnto.norb);
            nto.export_orbitals(*pr2, pnto.thresh);
            out << std::endl;
        }
        pr1->perform(density_type::particle, edm);
        pr1->perform(density_type::hole, hdm);
        if (m_h->is_active(wf_analysis_data_i::SA_NTO))
            add_to_average(edm, hdm);
    }
    else if (m_h->is_active(wf_analysis_data_i::NTO)) {

        orbital_params pnto = m_h->get_orbital_params(orbital_type::NTO);
        nto_analysis nto(s, c, tdm);
        nto.analyse(out, pnto.norb);
        nto.export_orbitals(*pr2, pnto.thresh);
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

            //store om_tot and om
            const ab_matrix &m_om = ct.omega();
            bool spin = true;
            if (m_om.is_alpha_eq_beta()) {
                spin = false;
            }
            frag_data_all[i][name].om_tot = ct.omega_total(spin);
            frag_data_all[i][name].om = m_om.alpha();

        }
    }

    if (m_h->is_active(wf_analysis_data_i::EXCITON)) {
        exciton_analysis(m_h->mom_builder(), tdm).analyse(out, 0);
    }
}


bool wf_analysis::setup_sa_ntos(std::ostream &out) {

    if (! m_h->is_active(wf_analysis_data_i::SA_NTO) || ! m_init_av)
        return false;
    if (m_sa.get())
        return true;

    ab_matrix edm(true), hdm(true);
    edm.alpha() = 0.5 * (m_edm_av.alpha() + m_edm_av.beta());
    hdm.alpha() = 0.5 * (m_hdm_av.alpha() + m_hdm_av.beta());

    const arma::mat &s = m_h->overlap();
    const ab_matrix &c = m_h->coefficients();

    nto_analysis sa_ntos(s, c, edm, hdm);
    if (m_h->is_active(wf_analysis_data_i::NTO)) {
        // Do not print out analysis of NTO spectrum here
        std::ostringstream dout;

        std::auto_ptr<orbital_printer_i> pr(m_h->orbital_printer("sa_nto",
                "State-averaged NTOs"));
        orbital_params pnto = m_h->get_orbital_params(orbital_type::NTO);
        sa_ntos.analyse(dout, pnto.norb);
        sa_ntos.export_orbitals(*pr, pnto.thresh);
    }
    m_sa = std::auto_ptr<sa_nto_analysis>(new sa_nto_analysis(s, sa_ntos));

    return true;
}

bool wf_analysis::post_process_optdm(std::ostream &out, const ab_matrix &tdm, const std::string &name, const double &ener) {

    if (m_h->n_ctnum_analyses() != 0) {

        for (size_t i = 0; i < m_h->n_ctnum_analyses(); i++) {

            // store state name
            frag_data_all[i][name].state_name = name;

            //store energy
            frag_data_all[i][name].dE_eV = ener;

            //computer descriptor and store
            const ctnum_analysis_i &ca = m_h->ctnum_analysis(i);
            const double &om_tot = frag_data_all[i][name].om_tot;
            const mat &om = frag_data_all[i][name].om;
            frag_data_all[i][name].descriptor = ca.compute_desc(om_tot, om);

        }
    }

    if (m_sa.get() == 0) return false;

    out << "Decomposition into state-averaged NTOs" << std::endl;
    m_sa->analyse(out, tdm);
    return true;
}


void wf_analysis::export_optdm(int prec) {

    if (m_h->n_ctnum_analyses() != 0 && !frag_data_all.empty()) {

        for (const auto& i : frag_data_all) {

            // file name
            std::string fname = "tden_summ_" + std::to_string(i.first + 1) + ".txt";

            // open file
            std::ofstream out;
            out.open(fname.c_str());

            // set precision
            out << std::setprecision(prec) << std::fixed;

            // header
            std::string header ("State          dE(eV)       f        ");
            out << header;

            // header: list of descriptor names
            auto itr = frag_data_all[i.first].begin();
            for (const auto& desc : itr->second.descriptor) {

                out << desc.first << std::setw(7);

            }
            out << std::endl;

            // dash line
            out << std::string(header.size() + itr->second.descriptor.size() * 10, '-') << std::endl;

            // data
            for (const auto& state : i.second) {

                out << state.second.state_name << std::setw(7);
                out << state.second.dE_eV << std::setw(2);
                out << state.second.f << std::setw(2);
                out << state.second.om_tot << std::setw(2);

                for (const auto& desc : state.second.descriptor) {

                    out << desc.second << std::setw(2);

                }
                out << std::endl;

            }

            out.close();

        }
    }
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


std::auto_ptr<wf_analysis> wf_analysis_static::analysis(0);

} // namespace libwfa
