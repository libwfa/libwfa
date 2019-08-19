#include <libwfa/libwfa.h>
#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

using namespace libwfa;

int main(int argc, char** argv)
{
#ifdef LIBWFA_DEBUG
    std::cout << "*** Debug mode activated ***" << std::endl << std::endl;
#endif

    char inp = ' ';
    char *pinp = &inp;
    molcas_wf_analysis_data *wfdata = libwfa::molcas_setup_wf_analysis_data(pinp);
    molcas_wf_analysis wf(wfdata);
    wf.run_analysis();

    return 0;
} // main
