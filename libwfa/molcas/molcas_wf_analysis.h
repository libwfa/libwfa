#ifndef LIBWFA_MOLCAS_WF_ANALYSIS_H
#define LIBWFA_MOLCAS_WF_ANALYSIS_H

#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

namespace libwfa {

class molcas_wf_analysis : public wf_analysis {
public:
    static const char k_clazz[]; //!< Class name
    
public:
    molcas_wf_analysis(wf_analysis_data_i *h, H5::H5File &file) : wf_analysis(h) { }
    
    void scf_analysis(molcas_wf_analysis_data *wfdata);
    void rasscf_analysis(molcas_wf_analysis_data *wfdata);
};
    
} // namespace libwfa

#endif // LIBWFA_MOLCAS_WF_ANALYSIS_H