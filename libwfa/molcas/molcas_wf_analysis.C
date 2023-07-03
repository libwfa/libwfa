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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <libwfa/libwfa.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include "H5Cpp.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis::k_clazz[] = "molcas_wf_analysis";

void molcas_wf_analysis::run_analysis() {

    std::string molcas_module = m_mdata->molcas_module();

    if (molcas_module=="SCF") {
        scf_analysis();
    }
    else if (molcas_module=="RASSCF") {
        rasscf_analysis(m_mdata->input()->refstate);
    }
    else if (molcas_module=="RASSI") {
        rassi_analysis(m_mdata->input()->refstate);
    }
    else {
        std::ostringstream os;
        os << std::endl << "Unsupported MOLCAS module: " << molcas_module;
        const std::string errmsg = os.str();
        throw libwfa_exception("main", "main", __FILE__, __LINE__, errmsg.c_str());
    }
}

void molcas_wf_analysis::scf_analysis() {
    header1("SCF MO Analysis");
    ab_matrix dm0 = m_mdata->build_dm(0, 0, true);
    analyse_opdm_ai("GS", "SCF ground state", dm0);
}

void molcas_wf_analysis::rasscf_analysis(size_t refstate) {
    header1("RASSCF Density Matrix Analysis");
    std::string label = m_mdata->rasscf_label();
    std::string labele = label;
    labele.erase(std::remove(labele.begin(), labele.end(), ' '), labele.end());

    // Density matrix
    arma::cube dens = m_mdata->read_cube_h5("DENSITY_MATRIX");
    if (refstate >= dens.n_slices)
        throw libwfa_exception(k_clazz, "rasscf_analysis", __FILE__, __LINE__, "refstate > nstate");

    double *dens_buf = dens.memptr();

    size_t dens_offs = dens.n_rows * dens.n_cols;

    // Spin-density matrix
    arma::cube sdens = m_mdata->read_cube_h5("SPINDENSITY_MATRIX");
    double *sdens_buf = sdens.memptr();

    double *smin = std::min_element(sdens_buf, sdens_buf + sdens.size());
    double *smax = std::max_element(sdens_buf, sdens_buf + sdens.size());
    bool aeqb_dens = true;
    if (*smax-*smin > 1.e-6) {
        std::cout << "Found non-vanishing spin-density, activating spin-analysis." << std::endl << std::endl;
        aeqb_dens = false;
    }

    ab_matrix dm0 = m_mdata->build_dm(dens_buf + refstate*dens_offs, sdens_buf + refstate*dens_offs, aeqb_dens);

    // Read the energies
    arma::vec ener = m_mdata->read_vec_h5("ROOT_ENERGIES");
    ener -= ener(refstate);

    // Loop over all states
    for (int istate = 0; istate < dens.n_slices; istate++) {
        std::ostringstream name, descr, header;
        descr << istate+1 << " " << label;
        name << istate+1 << labele;

        if (istate == refstate) {
            header << "RASSCF analysis for reference state " << descr.str();
            header2(header.str());
            m_mdata->energy_print(ener(istate), std::cout);

            analyse_opdm_ai(name.str(), descr.str(), dm0);
        }
        else {
            header << "RASSCF analysis for state " << descr.str();
            header2(header.str());
            m_mdata->energy_print(ener(istate), std::cout);

            ab_matrix ddm = m_mdata->build_dm(dens_buf + istate*dens_offs, sdens_buf + istate*dens_offs, aeqb_dens);
            ddm -= dm0;
            analyse_opdm_ai(name.str(), descr.str(), ddm, dm0);
        }
    }
}

void molcas_wf_analysis::rassi_analysis(size_t refstate) {
    header1("RASSI (Transition) Density Matrix Analysis");

    // Read the densities
    arma::cube tden = m_mdata->read_cube_h5("SFS_TRANSITION_DENSITIES");
    if (refstate >= tden.n_slices)
        throw libwfa_exception(k_clazz, "rassi_analysis", __FILE__, __LINE__, "refstate > nstate");

    double *tden_buf = tden.memptr();

    // Read the spin densities
    arma::cube tsden = m_mdata->read_cube_h5("SFS_TRANSITION_SPIN_DENSITIES");
    double *tsden_buf = tsden.memptr();

    //ab_matrix dm0 = m_mdata->build_dm_ao(tden_buf + (refstate + refstate * tden.n_cols)*tden.n_rows, tsden_buf + (refstate + refstate * tden.n_cols)*tden.n_rows, tden.n_rows);
    ab_matrix dm0 = m_mdata->build_dm_ao(tden_buf + (refstate + refstate * tden.n_cols)*tden.n_rows, NULL, tden.n_rows);

    // Read the energies
    arma::vec ener = m_mdata->read_vec_h5("SFS_ENERGIES");
    ener -= ener(refstate);

    // Read multiplicities and irreps; assign labels to the states
    int mult  [tden.n_slices];
    int irrep [tden.n_slices];
    std::vector<std::string> state_labels = m_mdata->rassi_labels(mult, irrep, tden.n_slices);

    // Read the transition moments
    arma::cube edip = m_mdata->read_cube_h5("SFS_EDIPMOM");

    // Loop over all transition densities
    for (int istate = 0; istate < tden.n_slices; istate++) {
        if (istate == refstate) {
            std::ostringstream descr, header;

            descr << state_labels[istate] << " " << std::setprecision(5) << ener(istate);

            header << "RASSI analysis for reference state " << state_labels[istate];
            header2(header.str());
            m_mdata->energy_print(ener(istate), std::cout);

            analyse_opdm_ai(state_labels[istate], descr.str(), dm0);
        }
        else {
            { // State/difference density analysis
                std::ostringstream descr, header;

                descr << state_labels[istate] << " " << std::setprecision(5) << ener(istate);

                header << "RASSI analysis for state " << state_labels[istate];
                header2(header.str());
                m_mdata->energy_print(ener(istate), std::cout);

                std::cout << "Using spin-traced density matrices ..." << std::endl << std::endl;

                const double *itden_buf  = tden_buf + (istate + istate * tden.n_cols)*tden.n_rows;

                // The state spin-densities are not understood consistently.
                //   Probably because of the contributions from inactive orbitals.
                //   They are ignored here
                const double *itsden_buf = NULL; //tsden_buf + (istate + istate * tden.n_cols)*tden.n_rows;
                //const double *itsden_buf = tsden_buf + (istate + istate * tden.n_cols)*tden.n_rows;

                ab_matrix ddm = m_mdata->build_dm_ao(itden_buf, itsden_buf, tden.n_rows);
                ddm -= dm0;

                analyse_opdm_ai(state_labels[istate], descr.str(), ddm, dm0);
            }
            { // Transition density analysis
                std::ostringstream name, descr, header;

                name << state_labels[refstate] << "-" << state_labels[istate];
                descr << name.str() << " " << std::setprecision(5) << ener(istate);

                header << "RASSI analysis for transiton from state " << refstate+1 << " to " << istate+1 << " (" << name.str() << ")";
                header2(header.str());
                m_mdata->energy_print(ener(istate), std::cout);

                int jstate = std::min(istate, (int)refstate);
                int kstate = std::max(istate, (int)refstate);

                const double *itden_buf  = tden_buf + (kstate + jstate * tden.n_cols)*tden.n_rows;
                const double *itsden_buf = tsden_buf + (kstate + jstate * tden.n_cols)*tden.n_rows;
                ab_matrix tdm = m_mdata->build_dm_ao(itden_buf, itsden_buf, tden.n_rows, irrep[refstate], irrep[istate]);
                if (istate > (int)refstate) // Transpose if the indices are switched
                    tdm.inplace_trans();

                const double energy = constants::au2eV * ener(istate);
                double osc = 0.;
                if (mult[kstate] == mult[jstate]) {
                    const double edip2 =
                        edip(jstate, kstate, 0) * edip(jstate, kstate, 0) +
                        edip(jstate, kstate, 1) * edip(jstate, kstate, 1) +
                        edip(jstate, kstate, 2) * edip(jstate, kstate, 2);
                    osc = 2./3. * ener(istate) * edip2;
                }

                analyse_optdm_ai(name.str(), descr.str(), tdm, energy, osc);
            }
        }
    }

    print_summary(std::cout);
    if (m_mdata->input()->add_info) add_molcas_info_fda();
    print_om_frag(std::cout);
}

void molcas_wf_analysis::header1(std::string title) {
    int lspace = (76 - title.size()) / 2;

    std::cout << "  " << std::string(76, '-') << std::endl;
    std::cout << std::string(lspace, ' ') << title << std::endl;
    std::cout << "  " << std::string(76, '-') << std::endl << std::endl;
}

void molcas_wf_analysis::header2(std::string title) {
    std::cout << std::endl << "  " << title << std::endl;
    std::cout << "  " << std::string(title.size(), '-') << std::endl;
}

void molcas_wf_analysis::analyse_opdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &ddm, const ab_matrix &dm0) {

    std::stringstream out;
    analyse_opdm(out, name, desc, ddm, dm0);
    std::cout << out.str();

    if (m_mdata->input()->add_info) add_molcas_info(out, "WFA1DDM");
}

void molcas_wf_analysis::analyse_opdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &sdm) {

    std::stringstream out;
    analyse_opdm(out, name, desc, sdm);
    std::cout << out.str();

    if (m_mdata->input()->add_info) add_molcas_info(out, "WFA1DM");
}

void molcas_wf_analysis::analyse_optdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &tdm, const double energy, const double osc) {

    std::stringstream out;
    analyse_optdm(out, name, desc, tdm, energy, osc);
    post_process_optdm(out, tdm);
    std::cout << out.str();

    if (m_mdata->input()->add_info) add_molcas_info(out, "WFA1TDM");
}

void molcas_wf_analysis::add_molcas_info(std::stringstream &in, std::string LAB) {
    size_t prec = 6;

    std::ofstream finfo;
    finfo.open("molcas_info", std::ofstream::app);

    double x;
    std::string str;

    finfo << std::setw(6);
    while (in >> str) {
        if (std::stringstream(str) >> x)
            if (x*x > 2.e-12)
            {
                finfo <<          LAB << "[" << m_info << "]=\"" << std::setprecision(prec)
                    << std::fixed << x << "\"" << std::endl;
                finfo << "#> " << LAB << "[" << m_info << "]=\"" << std::setprecision(prec)
                    << std::fixed << x << "\"/" << prec << std::endl;
                m_info += 1;
            }
    }
    finfo << "export " << LAB << std::endl;
}

void molcas_wf_analysis::add_molcas_info_fda() {
    size_t prec = 8;
    size_t i_info;

    std::ofstream finfo;
    finfo.open("molcas_info", std::ofstream::app);
    auto descs = m_mdata->prop_list();

    if (!frag_data_all.empty()) {
        for (const auto& desc : descs) {
            i_info = 0;
            for (const auto& i : frag_data_all) {
                for (const auto& state : i.second) {
                    auto descriptor = state.descriptor;
                    finfo <<          desc << "[" << i_info << "]=\"" << std::setprecision(prec)
                        << std::fixed << descriptor[desc] << "\"" << std::endl;
                    finfo << "#> " << desc << "[" << i_info << "]=\"" << std::setprecision(prec)
                        << std::fixed << descriptor[desc] << "\"/" << prec << std::endl;
                    i_info += 1;
                }
            }
            finfo << "export " << desc << std::endl;
        }
    }
}

} // namespace libwfa
