#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include "H5Cpp.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis::k_clazz[] = "molcas_wf_analysis";

void molcas_wf_analysis::scf_analysis() {
    header1("SCF MO Analysis");
    ab_matrix dm0 = m_mdata->build_dm(0, 0, true);
    analyse_opdm(std::cout, "GS", "SCF ground state", dm0);
}

void molcas_wf_analysis::rasscf_analysis(size_t refstate) {
    header1("RASSCF Density Matrix Analysis");

    // Density matrix
    arma::cube dens = m_mdata->read_dens_raw("DENSITY_MATRIX");
    double *dens_buf = dens.memptr();

    size_t dens_offs = dens.n_rows * dens.n_cols;

    // Spin-density matrix
    arma::cube sdens = m_mdata->read_dens_raw("SPINDENSITY_MATRIX");
    double *sdens_buf = sdens.memptr();

    double *smin = std::min_element(sdens_buf, sdens_buf + sdens.size());
    double *smax = std::max_element(sdens_buf, sdens_buf + sdens.size());
    bool aeqb_dens = true;
    if (*smax-*smin > 1.e-6) {
        std::cout << "Found non-vanishing spin-density, activating spin-analysis." << std::endl << std::endl;
        aeqb_dens = false;
    }

    ab_matrix dm0 = m_mdata->build_dm(dens_buf + refstate*dens_offs, sdens_buf + refstate*dens_offs, aeqb_dens);

    // Loop over all states
    for (int istate = 0; istate < dens.n_slices; istate++) {
        if (istate == refstate) {
            std::ostringstream name, header;
            name << "S_" << refstate+1;
            header << "RASSCF analysis for reference state " << name.str();
            header2(header.str());

            analyse_opdm(std::cout, name.str(), name.str(), dm0);
        }
        else {
            std::ostringstream name, header;
            name << "S_" << istate+1;
            header << "RASSCF analysis for state " << name.str();
            header2(header.str());
    
            ab_matrix ddm = m_mdata->build_dm(dens_buf + istate*dens_offs, sdens_buf + istate*dens_offs, aeqb_dens);
            ddm -= dm0;
            analyse_opdm(std::cout, name.str(), name.str(), ddm, dm0);
        }
    }
}

void molcas_wf_analysis::rassi_analysis(size_t refstate) {
    header1("RASSI Transition Density Matrix Analysis");

    arma::cube tden = m_mdata->read_dens_raw("SFS_TRANSITION_DENSITIES");
    const double *tden_buf = tden.memptr();
    
    ab_matrix dm0 = m_mdata->build_dm_ao(tden_buf + (refstate + refstate * tden.n_cols)*tden.n_rows, tden.n_rows);
    
    // Loop over all transition densities
    for (int istate = 0; istate < tden.n_slices; istate++) {
        if (istate == refstate) {
            std::ostringstream name, header;
            name << "S_" << refstate+1;
            header << "RASSI analysis for reference state " << std::setw(3) << refstate+1;
            header2(header.str());
            analyse_opdm(std::cout, name.str(), name.str(), dm0);
        }
        else {
            { // State/difference density analysis
                std::ostringstream name, header;
                name << "S_" << istate+1;
                header << "RASSI analysis for state " << std::setw(3) << istate+1;
                header2(header.str());
                
                const double *itden_buf = tden_buf + (istate + istate * tden.n_cols)*tden.n_rows;
                ab_matrix ddm = m_mdata->build_dm_ao(itden_buf, tden.n_rows);
                ddm -= dm0;
                analyse_opdm(std::cout, name.str(), name.str(), ddm, dm0);
            }
            { // Transition density analysis
                int jstate = std::min(istate, (int)refstate);
                int kstate = std::max(istate, (int)refstate);
                std::ostringstream name, header;                
                name << "T_" << jstate + 1 << "-" << kstate + 1;
                header << "RASSI analysis for transiton from " << jstate+1 << " to " << kstate+1;
                header2(header.str());
                
                const double *itden_buf = tden_buf + (kstate + jstate * tden.n_cols)*tden.n_rows;
                ab_matrix tdm = m_mdata->build_dm_ao(itden_buf, tden.n_rows);
                tdm.inplace_trans();
                analyse_optdm(std::cout, name.str(), name.str(), tdm);
            }
        }
    }
}

void molcas_wf_analysis::header1(std::string title) {
    int lspace = (76 - title.size()) / 2;

    std::cout << "  " << std::string(76, '-') << std::endl;
    std::cout << std::string(lspace, ' ') << title << std::endl;
    std::cout << "  " << std::string(76, '-') << std::endl << std::endl;
}

void molcas_wf_analysis::header2(std::string title) {
    std::cout << "  " << title << std::endl;
    std::cout << "  " << std::string(title.size(), '-') << std::endl;
}

}