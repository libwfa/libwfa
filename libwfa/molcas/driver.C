#include <libwfa/libwfa.h>
#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

using namespace libwfa;

int wfa_driver(char *inp)
{
#ifdef LIBWFA_DEBUG
    std::cout << "*** Debug mode activated ***" << std::endl << std::endl;
#endif

    molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(inp);
    molcas_wf_analysis wf(wfdata);
    wf.run_analysis();

    return 0;
} // main

extern "C"
{
    void wfa_driver_(int *rc, char *inp, int linp);
}

void wfa_driver_(int *rc, char *inp, int linp)
{
    inp[linp--] = '\0'; // NULL terminate the string
    *rc = wfa_driver(inp);
}
