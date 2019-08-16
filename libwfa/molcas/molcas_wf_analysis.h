#ifndef LIBWFA_MOLCAS_WF_ANALYSIS_H
#define LIBWFA_MOLCAS_WF_ANALYSIS_H

#include <libwfa/wf_analysis.h>
#include <libwfa/molcas/molcas_wf_analysis_data.h>
#include "libwfa/core/constants.h"

namespace libwfa {

class molcas_wf_analysis : public wf_analysis {
private:
    molcas_wf_analysis_data *m_mdata;
    size_t m_info; //!< Index for add_molcas_info

public:
    static const char k_clazz[]; //!< Class name

public:
    molcas_wf_analysis(molcas_wf_analysis_data *h) : wf_analysis(h), m_mdata(h), m_info(0) { }

    void run_analysis();

private:
    void scf_analysis();
    void rasscf_analysis(size_t refstate);
    void rassi_analysis(size_t refstate);
    void header1(std::string title);
    void header2(std::string title);

    /** \brief Analysis of sdm and ddm and add_info
     **/
    void analyse_opdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &ddm, const ab_matrix &dm0);

    /** \brief Analysis of sdm and add_info
     **/
    void analyse_opdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &sdm);

    /** \brief Analysis of tdm and add_info
     **/
    void analyse_optdm_ai(const std::string &name, const std::string &desc,
        const ab_matrix &tdm, const double energy, const double osc);

    /** \brief Append to molcas_info file

        Add information to the molcas_info file in the form expected
        by the molcas verify code.
     **/
    void add_molcas_info(std::stringstream &out);
};

} // namespace libwfa

#endif // LIBWFA_MOLCAS_WF_ANALYSIS_H
