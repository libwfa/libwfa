#include "H5Cpp.h"
#include <libwfa/molcas/molcas_wf_analysis.h>
#include <libwfa/libwfa.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

using namespace H5;
using namespace libwfa;

int main(int argc, char** argv)
{
    H5std_string file_name("molcas.rasscf.h5");
    H5std_string *ref_name = &file_name;
    int ref_state = 1;
    if (argc>=2){
        file_name = argv[1];
        if (argc>=3) {
            ref_state = atoi(argv[2]);
            if (argc>=4) {
                ref_name = new H5std_string(argv[3]);
            }
        }
    }
    std::cout << "Starting analysis of Molcas HDF5 file " << file_name << std::endl;

    H5File file( file_name, H5F_ACC_RDWR );

    molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(file);
    molcas_wf_analysis wf(wfdata);

    // Check what kind of job was performed
    H5std_string molcas_module;
    {
        Attribute Att = file.openGroup("/").openAttribute("MOLCAS_MODULE");
        StrType strtype(PredType::C_S1, 16);
        Att.read(strtype, molcas_module);
    }

    if (molcas_module=="SCF") {
        wf.header1("SCF MO Analysis");
        wf.scf_analysis();
    }
    else if (molcas_module=="RASSCF") {
        wf.header1("RASSCF Density Matrix Analysis");
        std::cout << "Using state " << ref_state << " of file " <<
            *ref_name << " as a reference for attachment/detachment analysis." << std::endl;
        wf.rasscf_analysis();
    }
    else if (molcas_module=="RASSI") {
        wf.header1("RASSI Transition Density Matrix Analysis");
        wf.rassi_analysis();
    }
    else {
        std::ostringstream os;
        os << std::endl << "Unknown MOLCAS module: " << molcas_module;
        const std::string errmsg = os.str();
        throw libwfa_exception("main", "main", __FILE__, __LINE__, errmsg.c_str());
    }
    
    return 0;
} // main