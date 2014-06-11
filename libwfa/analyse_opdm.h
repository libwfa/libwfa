#ifndef LIBWFA_ANALYSE_OPDM_H
#define LIBWFA_ANALYSE_OPDM_H

#include <map>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>
#include <libwfa/export/pop_printer_i.h>

namespace libwfa {

/** \brief Combines various analyses of excited state density matrices

    \ingroup libwfa
 **/
class analyse_opdm {
public:
    /** \brief Flag when to perform specific population analysis
     **/
    typedef enum {
        pa_none = 0, //!< Do not perform PA
        pa_dm = 1, //!< Perform PA on state density matrix
        pa_ad = 2, //!< Perform PA on attachment / detachment density matrices
        pa_both = 3 //!< Perform both PAs
    } pa_flag;

private:
    struct pa {
        const pop_analysis_i &analysis;
        const pop_printer_i &printer;
        pa_flag flag;

        pa(const pop_analysis_i &a, const pop_printer_i &p, pa_flag fl) :
            analysis(a), printer(p), flag(fl) { }
    };
    typedef std::map<std::string, pa> pa_map_t;

private:
    pa_map_t m_pa; //!< List of population analyses
    const ev_printer_i *m_pr[2]; //!< Formating objects for NO and NDO summary

    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficient matrix

    const ab_matrix &m_dm1; //!< State or difference density matrix
    std::auto_ptr<ab_matrix> m_dm2; //!< State or differnce density matrix
    const ab_matrix &m_sdm; //!< State density matrix
    const ab_matrix &m_ddm; //!< Difference density matrix

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param dm State density matrix
     **/
    analyse_opdm(const arma::Mat<double> &s,
        const ab_matrix &c, const ab_matrix &dm);

    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param dm0 Ground state density matrix
        \param dm State or difference density matrix
        \param is_diff Is difference density?
     **/
    analyse_opdm(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &dm0, const ab_matrix &dm, bool is_diff = true);

    /** \brief Register orbital analysis and printer
        \param pr Printer for orbital analysis
        \param is_no True: register NO analysis; false: register NDO analysis
     **/
    void do_register(const ev_printer_i &pr, bool is_no);

    /** \brief Register population analyses that should be performed
        \param name Name for population analysis
        \param pa Population analysis
        \param pr Population printer / formatter
        \param fl Flag for which density the population analysis is used
     **/
    void do_register(const std::string &name, const pop_analysis_i &pa,
        const pop_printer_i &pr, pa_flag fl);

    /** \brief Performs density matrix analyses
        \param pr Density and orbital export / printer
        \param out Output stream

        Perform the following analyses:
        - NO analysis
        - NDO analysis
        - Population analysis
     **/
    void perform(export_data_i &pr, std::ostream &out) const;

private:
    static std::auto_ptr<ab_matrix> build_dm(const ab_matrix &dm,
        const ab_matrix &dm0, bool is_diff);
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_OPDM_H
