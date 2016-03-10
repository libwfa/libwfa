#ifndef LIBWFA_MOLCAS_WF_ANALYSIS_H
#define LIBWFA_MOLCAS_WF_ANALYSIS_H

#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>

namespace libwfa {

class molcas_wf_analysis : public wf_analysis {
private:
    molcas_wf_analysis_data *m_mdata;
    
public:
    static const char k_clazz[]; //!< Class name
    
public:
    molcas_wf_analysis(molcas_wf_analysis_data *h) : wf_analysis(h), m_mdata(h) { }
    
    void scf_analysis();
    void rasscf_analysis(size_t refstate);
    void rassi_analysis(size_t refstate);

private:
    void header1(std::string title);
    void header2(std::string title);
};
    
} // namespace libwfa

#endif // LIBWFA_MOLCAS_WF_ANALYSIS_H