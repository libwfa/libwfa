#ifndef LIBWFA_WF_ANALYSIS_H
#define LIBWFA_WF_ANALYSIS_H

#include <libwfa/analyses/sa_nto_analysis.h>
#include <set>
#include "wf_analysis_data_i.h"

namespace libwfa {


/** \brief General parameters for analyses

    Stores the orbital parameters for NO, NDO, and NTO analyses.
    The orbital parameters consist of the number of leading orbital pairs
    (integer) to print and a threshold for important orbitals to export.
 **/
struct wfa_params {
    std::pair<size_t, double> no, ndo, nto; //!< Orbital parameters

    wfa_params() : no(4, 1e-3), ndo(5, 1e-3), nto(5, 1e-3) { }
};

/** \brief Wave function analysis class

    Wrapper class to perform all available analyses.

    If new analyses are added to the library, please amend the respective
    functions below.

    Currently the following analyses can be activated:
    - NO analysis: no
    - NDO analysis: ndo
    - Form attachment/detachment densities: ad
    - Exciton analysis based on a/d densities: exciton_ad
    - NTO analysis: nto
    - Formation of electron/hole densities: eh
    - Exciton analysis based on transition densities: exciton
    - State-averaged NTO analysis: sa_nto

    Additionally, population and CT number analyses will be performed, if
    provided by the analysis data object.

    \ingroup libwfa
 **/
class wf_analysis {
private:
    typedef enum {
        NO = 0,     //!< NO analysis
        NDO,        //!< NDO analysis
        FORM_AD,    //!< Form attachment/detachment densities
        EXCITON_AD, //!< Exciton analysis on a/d densities
        NTO,        //!< NTO analysis
        FORM_EH,    //!< Form electron/hole densities
        EXCITON,    //!< Exciton analysis on transition density
        SA_NTO,     //!< State-averaged NTO analysis
        NA          //!< Number of analyses
    } ana_type;

private:
    std::auto_ptr<wf_analysis_data_i> m_h; //!< Analysis data
    std::bitset<NA> m_p1; //!< Analyses to run
    wfa_params m_p2;

    std::auto_ptr<sa_nto_analysis> m_sa; //!< State-averaged NTO analysis
    ab_matrix m_edm_av; //!< Averaged electron density
    ab_matrix m_hdm_av; //!< Averaged hole density
    bool m_init_av; //!< Whether the above are initialized?

public:
    /** \brief Initializes the wave function analysis
        \param h Analysis data
        \param p Parameters for analyses
     **/
    wf_analysis(wf_analysis_data_i *h, const wfa_params &p = wfa_params()) :
        m_h(h), m_p2(p), m_sa(0), m_init_av(false) {
        m_p1.set();
    }

    /** \brief Initializes the wave function analysis
        \param h Analysis data
        \param p1 Analyses to activate
        \param p2 Parameters for analyses
     **/
    wf_analysis(wf_analysis_data_i *h, const std::set<std::string> &p1,
            const wfa_params &p2 = wfa_params());

    /** \brief Activate a specific analysis
     **/
    void activate(const std::string &a);

    /** \brief Deactivate a specific analysis
     **/
    void deactivate(const std::string &a);

    /** \brief Test if a specific analysis is active
     **/
    bool is_active(const std::string &a) const {
        ana_type p = convert(a);
        return (p < NA ? m_p1.test(p) : false);
    }

    /** \brief Perform analysis of state and difference density matrix
        \param out Output stream
        \param name Name of state (useable as filename)
        \param desc Description of state (one-line comment)
        \param ddm Difference density matrix in AO
        \param dm0 Ground state density matrix in AO
     **/
    void analyse_opdm(std::ostream &out,
        const std::string &name, const std::string &desc,
        const ab_matrix &ddm, const ab_matrix &dm0);

    /** \brief Perform analysis of state density matrix
            \param out Output stream
            \param name Name of state (useable as filename)
            \param desc Description of state (one-line comment)
            \param sdm State density matrix in AO
     **/
    void analyse_opdm(std::ostream &out, const std::string &name,
        const std::string &desc, const ab_matrix &sdm);

    /** \brief Perform analysis of transition density matrix
        \param out Output stream
        \param name Name of state (useable as filename)
        \param desc Description of state (one-line comment)
        \param tdm Transition density matrix in AO
     **/
    void analyse_optdm(std::ostream &out, const std::string &name,
        const std::string &desc, const ab_matrix &tdm);

    /** \brief Constructs state-averaged NTOs and sets up the analysis
        \param out Output stream
        \return True, if successful
     **/
    bool setup_sa_ntos(std::ostream &out);

    /** \brief Perform SA-NTO analysis of transition density matrix
        \param out Output stream
        \param tdm Transition density matrix in AO
        \return True, if successful
     **/
    bool post_process_optdm(std::ostream &out, const ab_matrix &tdm);

    /** \brief Reset to original state (at construction)
     **/
    void reset() {
        m_init_av = false;
        delete m_sa.release();
    }

private:
    void add_to_average(const ab_matrix &edm, const ab_matrix &hdm);
    static ana_type convert(const std::string &a);
};


/** \brief Holder for static analysis object

    The analysis object needs to be initialized, for the C-style interface to
    work properly.
 **/
struct wf_analysis_static {
    static std::auto_ptr<wf_analysis> analysis; //!< Static analysis object
};


} // namespace libwfa

#endif // LIBWFA_WF_ANALYSIS_H
