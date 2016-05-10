// #include "H5Cpp.h"
#include <libwfa/libwfa.h>
#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

// using namespace H5;
using namespace libwfa;

int wfa_driver()
{
#ifdef LIBWFA_DEBUG
    std::cout << "*** Debug mode activated ***" << std::endl << std::endl;
#endif

    molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data();
    molcas_wf_analysis wf(wfdata);
    wf.run_analysis();

    return 0;
} // main

extern "C"
{
    void wfa_driver_(int *rc);
}

void wfa_driver_(int *rc)
{
    *rc = wfa_driver();
}