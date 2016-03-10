#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include "H5Cpp.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis::k_clazz[] = "molcas_wf_analysis";

void molcas_wf_analysis::scf_analysis() {

    ab_matrix dm0 = m_mdata->build_dm(0, 0, true);
    analyse_opdm(std::cout, "GS", "SCF ground state", dm0);
}

void molcas_wf_analysis::rasscf_analysis() {

    // Density matrix
    arma::cube dens = m_mdata->read_dens_raw("DENSITY_MATRIX");
    double *dens_buf = dens.memptr();

    size_t nexc = dens.n_slices - 1;
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

    header2("RASSCF ground state");
    ab_matrix dm0 = m_mdata->build_dm(dens_buf, sdens_buf, aeqb_dens);
    analyse_opdm(std::cout, "GS", "RASSCF ground state", dm0);

    // Loop over all states
    for (int istate = 1; istate <= nexc; istate++) {
        std::ostringstream name, descr;
        name << "ES_" << istate;
        descr << "RASSCF excited state " << std::setw(3) << istate;
        header2(descr.str());

        ab_matrix ddm = m_mdata->build_dm(dens_buf + istate*dens_offs, sdens_buf + istate*dens_offs, aeqb_dens);
        ddm -= dm0;
        analyse_opdm(std::cout, name.str(), descr.str(), ddm, dm0);
    }
}

void molcas_wf_analysis::rassi_analysis() {
    // Density matrix
    arma::cube tden = m_mdata->read_dens_raw("SFS_TRANSITION_DENSITIES");
    const double *tden_buf = tden.memptr();
    
    // Loop over all transition densities
    for (int iden = 0; iden < tden.n_slices; iden++) {
        for (int jden = 0; jden < tden.n_cols; jden++) {
            if (iden==jden) {
                std::ostringstream name, descr;
                name << "ES_" << iden;
                descr << "RASSI analysis for state " << std::setw(3) << iden+1;
                header2(descr.str());
                ab_matrix dm = m_mdata->build_dm_ao(tden_buf, tden.n_rows);
                dm *= 0.5;
                analyse_opdm(std::cout, name.str(), descr.str(), dm);
            }
            else if (jden > iden) {
                std::ostringstream name, descr;
                name << "T_" << iden + 1 << "-" << jden + 1;
                descr << "RASSI analysis for transiton from " << iden+1 << " to " << jden+1;
                header2(descr.str());
                
                ab_matrix tdm = m_mdata->build_dm_ao(tden_buf, tden.n_rows);
                tdm *= 0.5;
                tdm.inplace_trans();
                analyse_optdm(std::cout, name.str(), name.str(), tdm);
            }
            
            tden_buf += tden.n_rows;
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