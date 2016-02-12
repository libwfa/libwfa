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
    
    H5File file( file_name, H5F_ACC_RDONLY ); // read-only for now
    
    libwfa::molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(file);
    libwfa::wf_analysis wf(wfdata);
    
    // First do density matrix analysis
    std::cout << "  " << std::string(76, '-') << std::endl;
    std::cout << std::string(28, ' ') << "Density Matrix Analysis" << std::endl;
    std::cout << "  " << std::string(76, '-') << std::endl << std::endl;
    {
        DataSet Set = file.openDataSet("DENSITY_MATRIX");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 3)
            throw libwfa_exception("main", "main", __FILE__, __LINE__, "Inconsistent rank for DM");
        
        hsize_t dims[rank];
        Space.getSimpleExtentDims(dims, NULL);
        size_t nexc = dims[0]-1;        
        size_t dens_offs = dims[1] * dims[2];
        
        double dens_buf[dims[0] * dims[1] * dims[2]];
        Set.read(&dens_buf, PredType::NATIVE_DOUBLE);
        
/*        { // closed shell
            std::cout << "Closed shell" << std::endl;
            ab_matrix dm0(44, 44);
            dm0.alpha().zeros();
            int nocc = 16;
            for (int iocc=0; iocc<nocc; iocc++)
                dm0.alpha().at(iocc, iocc) = 1.;
            
            arma::mat ca = wfdata->coefficients().alpha();
            dm0.alpha() = ca * dm0.alpha() * ca.t();
            
            std::ostringstream name;
            name << "hf";
            
            wf.analyse_opdm(std::cout, name.str(), name.str(), dm0);       
        }*/
        
        std::cout << "  Ground state:" << std::endl;
        std::cout << "  " << std::string(18, '-') << std::endl;        
        ab_matrix dm0 = wfdata->build_dm(dens_buf);
        wf.analyse_opdm(std::cout, "gs", "gs", dm0);
        
        // Loop over all states
        for (int istate = 1; istate <= nexc; istate++) {
            std::cout << "  Excited state " << std::setw(3) << istate
                    << ":" << std::endl;
            std::cout << "  " << std::string(18, '-') << std::endl;
            
            ab_matrix ddm = wfdata->build_dm(dens_buf + istate*dens_offs); ddm -= dm0;
            wf.analyse_opdm(std::cout, "es", "es", ddm, dm0);
        }
    } // density matrix analysis
}