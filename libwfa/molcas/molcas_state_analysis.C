#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "H5Cpp.h"
#include <libwfa/libwfa.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

using namespace H5;
using namespace libwfa;

int main(int argc, char** argv)
{
    H5std_string file_name("molcas.rasscf.h5");
    if (argc>=2){
        file_name = argv[1];
    }
    std::cout << "Starting analysis of Molcas HDF5 file " << file_name << std::endl;

    H5File file( file_name, H5F_ACC_RDWR ); // read-only for now

    libwfa::molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(file);
    libwfa::wf_analysis wf(wfdata);

    // Check what kind of job was performed
    H5std_string molcas_module;
    {
        Attribute Att = file.openGroup("/").openAttribute("MOLCAS_MODULE");
        StrType strtype(PredType::C_S1, 16);
        Att.read(strtype, molcas_module);
    }

    if (molcas_module=="SCF") {
        std::cout << "  " << std::string(76, '-') << std::endl;
        std::cout << std::string(31, ' ') << "SCF MO Analysis" << std::endl;
        std::cout << "  " << std::string(76, '-') << std::endl << std::endl;

        ab_matrix dm0 = wfdata->build_dm(0, 0, true);
        wf.analyse_opdm(std::cout, "GS", "gs", dm0);
    }
    else if (molcas_module=="RASSCF") {
        std::cout << "  " << std::string(76, '-') << std::endl;
        std::cout << std::string(23, ' ') << "RASSCF Density Matrix Analysis" << std::endl;
        std::cout << "  " << std::string(76, '-') << std::endl << std::endl;

        // Density matrix
        DataSet Set = file.openDataSet("DENSITY_MATRIX");
        hsize_t dims[3];
        {
            DataSpace Space = Set.getSpace();
            if (Space.getSimpleExtentNdims() != 3)
                throw libwfa_exception("main", "main", __FILE__, __LINE__, "Inconsistent rank for DM");

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
            DataSet Set_s = file.openDataSet("SPINDENSITY_MATRIX");
            if (Set_s.getSpace().getSimpleExtentNdims() != 3)
                throw libwfa_exception("main", "main", __FILE__, __LINE__, "Inconsistent rank for Spin-DM");

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
        wf.analyse_opdm(std::cout, "GS", "gs++", dm0);

        // Loop over all states
        for (int istate = 1; istate <= nexc; istate++) {
            std::ostringstream name, descr;
            name << "ES_" << istate;
            descr << "Excited state " << std::setw(3) << istate;

            std::cout << "  " << descr.str() << ":" << std::endl;
            std::cout << "  " << std::string(18, '-') << std::endl;

            ab_matrix ddm = wfdata->build_dm(dens_buf + istate*dens_offs, sdens_buf + istate*dens_offs, aeqb_dens); ddm -= dm0;
            wf.analyse_opdm(std::cout, name.str(), descr.str(), ddm, dm0);
        }
    } // RASSCF
    else {
        std::ostringstream os;
        os << std::endl << "Unknown molcas MOLCAS module: " << molcas_module;
        const std::string errmsg = os.str();
        throw libwfa_exception("main", "main", __FILE__, __LINE__, errmsg.c_str());
    }
    return 0;
} // main