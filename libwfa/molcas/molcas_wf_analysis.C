#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include "H5Cpp.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis::k_clazz[] = "molcas_wf_analysis";

void molcas_wf_analysis::scf_analysis(molcas_wf_analysis_data *wfdata) {
    
    std::cout << "  " << std::string(76, '-') << std::endl;
    std::cout << std::string(31, ' ') << "SCF MO Analysis" << std::endl;
    std::cout << "  " << std::string(76, '-') << std::endl << std::endl;

    ab_matrix dm0 = wfdata->build_dm(0, 0, true);
    analyse_opdm(std::cout, "GS", "gs", dm0);
}

void molcas_wf_analysis::rasscf_analysis(molcas_wf_analysis_data *wfdata) {
    
    std::cout << "  " << std::string(76, '-') << std::endl;
    std::cout << std::string(23, ' ') << "RASSCF Density Matrix Analysis" << std::endl;
    std::cout << "  " << std::string(76, '-') << std::endl << std::endl;

    // Density matrix
    DataSet Set = wfdata->h5file().openDataSet("DENSITY_MATRIX");
    hsize_t dims[3];
    {
        DataSpace Space = Set.getSpace();
        if (Space.getSimpleExtentNdims() != 3)
            throw libwfa_exception(k_clazz, "rasscf_analysis", __FILE__, __LINE__, "Inconsistent rank for DM");

        Space.getSimpleExtentDims(dims, NULL);
    }
    size_t nexc = dims[0]-1;
    size_t dens_offs = dims[1] * dims[2];

    double dens_buf[dims[0] * dims[1] * dims[2]];
    Set.read(&dens_buf, PredType::NATIVE_DOUBLE);

    // Spin-density matrix
    bool aeqb_dens = true;
    double sdens_buf[dims[0] * dims[1] * dims[2]];
    {
        DataSet Set_s = wfdata->h5file().openDataSet("SPINDENSITY_MATRIX");
        if (Set_s.getSpace().getSimpleExtentNdims() != 3)
            throw libwfa_exception(k_clazz, "rasscf_analysis", __FILE__, __LINE__, "Inconsistent rank for Spin-DM");

        Set_s.read(&sdens_buf, PredType::NATIVE_DOUBLE);

        double *smin = std::min_element(sdens_buf, sdens_buf + dims[0] * dims[1] * dims[2]);
        double *smax = std::max_element(sdens_buf, sdens_buf + dims[0] * dims[1] * dims[2]);
        if (*smax-*smin > 1.e-6) {
            std::cout << "Found non-vanishing spin-density, activating spin-analysis." << std::endl << std::endl;
            aeqb_dens = false;
        }
    }

    std::cout << "  Ground state:" << std::endl;
    std::cout << "  " << std::string(18, '-') << std::endl;
    ab_matrix dm0 = wfdata->build_dm(dens_buf, sdens_buf, aeqb_dens);
    analyse_opdm(std::cout, "GS", "gs++", dm0);

    // Loop over all states
    for (int istate = 1; istate <= nexc; istate++) {
        std::ostringstream name, descr;
        name << "ES_" << istate;
        descr << "Excited state " << std::setw(3) << istate;

        std::cout << "  " << descr.str() << ":" << std::endl;
        std::cout << "  " << std::string(18, '-') << std::endl;

        ab_matrix ddm = wfdata->build_dm(dens_buf + istate*dens_offs, sdens_buf + istate*dens_offs, aeqb_dens); ddm -= dm0;
        analyse_opdm(std::cout, name.str(), descr.str(), ddm, dm0);
    }
}

}