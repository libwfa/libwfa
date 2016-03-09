#include "H5Cpp.h"
#include <libwfa/molcas/molcas_wf_analysis.h>
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

    molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(file);
    molcas_wf_analysis wf(wfdata, file);

    // Check what kind of job was performed
    H5std_string molcas_module;
    {
        Attribute Att = file.openGroup("/").openAttribute("MOLCAS_MODULE");
        StrType strtype(PredType::C_S1, 16);
        Att.read(strtype, molcas_module);
    }

    if (molcas_module=="SCF") {
        wf.scf_analysis(wfdata);
    }
    else if (molcas_module=="RASSCF") {
        wf.rasscf_analysis(wfdata);
    }
    else {
        std::ostringstream os;
        os << std::endl << "Unknown molcas MOLCAS module: " << molcas_module;
        const std::string errmsg = os.str();
        throw libwfa_exception("main", "main", __FILE__, __LINE__, errmsg.c_str());
    }
    return 0;
} // main