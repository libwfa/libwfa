//************************************************************************
//* This file is part of libwfa.                                         *
//*                                                                      *
//* libwfa is free software; you can redistribute and/or modify          *
//* it under the terms of the BSD 3-Clause license.                      *
//* libwfa is distributed in the hope that it will be useful, but it     *
//* is provided "as is" and without any express or implied warranties.   *
//* For more details see the full text of the license in the file        *
//* LICENSE.                                                             *
//*                                                                      *
//* Copyright (c) 2014, F. Plasser and M. Wormit. All rights reserved.   *
//* Modifications copyright (C) 2019, Loughborough University.           *
//************************************************************************


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
    std::unique_ptr<density_printer_i> pr1(m_h->density_printer(name, desc));
    std::unique_ptr<orbital_printer_i> pr2(m_h->orbital_printer(name, desc));

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
    std::unique_ptr<density_printer_i> pr1(h.density_printer(name, desc));
    std::unique_ptr<orbital_printer_i> pr2(h.orbital_printer(name, desc));

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
    const std::string &desc, const ab_matrix &tdm, double energy) {

    // Create printer for orbitals and densities
    std::unique_ptr<density_printer_i> pr1(m_h->density_printer(name, desc));
    std::unique_ptr<orbital_printer_i> pr2(m_h->orbital_printer(name, desc));

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
            std::unique_ptr<ctnum_printer_i> cpr(m_h->ctnum_printer(i, name, desc));

            out << cname << std::endl;
            ctnumbers ct(ca, tdm);
            ct.analyse(out);
            ct.do_export(*cpr);
            out << std::endl;

            //Compute spin-traced Omega matrices
            const ab_matrix &m_om = ct.omega();
            frag_data_all[i][name].om_tot = ct.omega_total(false) + ct.omega_total(true);
            const mat om_at = m_om.alpha() + m_om.beta();

            // Transform to fragments, compute descriptors and store
            frag_data_all[i][name].om_frag = ca.compute_omFrag(om_at);
            frag_data_all[i][name].descriptor = ca.compute_descriptors(frag_data_all[i][name].om_tot,
                frag_data_all[i][name].om_frag);

            // store state name
            frag_data_all[i][name].state_name = name;

            //store energy
            frag_data_all[i][name].dE_eV = energy;
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

        std::unique_ptr<orbital_printer_i> pr(m_h->orbital_printer("sa_nto",
                "State-averaged NTOs"));
        orbital_params pnto = m_h->get_orbital_params(orbital_type::NTO);
        sa_ntos.analyse(dout, pnto.norb);
        sa_ntos.export_orbitals(*pr, pnto.thresh);
    }
    m_sa = std::unique_ptr<sa_nto_analysis>(new sa_nto_analysis(s, sa_ntos));

    return true;
}

bool wf_analysis::post_process_optdm(std::ostream &out, const ab_matrix &tdm) {

    if (m_sa.get() == 0) return false;

    out << "Decomposition into state-averaged NTOs" << std::endl;
    m_sa->analyse(out, tdm);
    return true;
}


void wf_analysis::print_summary(std::ostream &out, const int &prec, const int &width) {

    if (m_h->n_ctnum_analyses() != 0 && !frag_data_all.empty()) {
        int hwidth;
        for (const auto& i : frag_data_all) {

            // set precision
            out << std::setprecision(prec) << std::fixed;

            // header with list of descriptor names
            std::string header ("State            dE(eV)      f         ");
            auto descs = m_h->prop_list();
            hwidth = header.size() + descs.size() * width;

            out  << std::endl;
            out << std::string(5, '=') << " TheoDORE-style Summary ";
            out  << std::string(hwidth - 27, '=') << std::endl;
            out << "| " << header;
            for (const auto& desc : descs) {
                out << std::left << std::setw(width) << desc;
            }
            out << std::endl;

            // dash line
            out << "| " << std::string(hwidth, '-') << std::endl;

            // data
            for (const auto& state : i.second) {
                out << "| ";
                out << std::left << std::setw(14)     << state.second.state_name;
                out << std::right << std::setw(width) << state.second.dE_eV;
                out << std::right << std::setw(width) << state.second.f;

                auto descriptors = state.second.descriptor;
                for (const auto& desc : descs) {
                    out << std::right << std::setw(width) << descriptors[desc];
                }
                out << std::endl;
            }
        }
        out << std::string(hwidth + 2, '=') << std::endl << std::endl;
    }
}

void wf_analysis::print_om_frag(std::ostream &out, const std::string ofile) {
    if (m_h->n_ctnum_analyses() == 0 || frag_data_all.empty())
        return;

    std::cout << " Writing fragment Omega matrix to " << ofile << std::endl;
    std::ofstream fout;
    fout.open(ofile.c_str());

    //fout << at_lists.size();
    bool header = false;
    for (const auto& i : frag_data_all) {
        for (const auto& state : i.second) {
            int nfrag = size(state.second.om_frag)[0];

            if (!header) {
                fout << nfrag << std::endl;
                header = true;
            }

            fout << state.second.state_name << " ";
            fout << " " << state.second.om_tot;

            for (size_t ifrag = 0; ifrag < nfrag; ifrag++) {
                for (size_t jfrag = 0; jfrag < nfrag; jfrag++) {
                    fout << " " << state.second.om_frag(ifrag, jfrag);
                }
            }

            fout << std::endl;
        }
    }

    fout.close();
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


std::unique_ptr<wf_analysis> wf_analysis_static::analysis(nullptr);

} // namespace libwfa
