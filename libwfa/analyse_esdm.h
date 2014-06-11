#ifndef LIBWFA_ANALYSE_ESDM_H
#define LIBWFA_ANALYSE_ESDM_H

#include <map>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/export/pop_printer_i.h>

namespace libwfa {

/** \brief Combines various analyses of excited state density matrices

    \ingroup libwfa
 **/
class analyse_esdm {
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
    pa_map_t m_lst; //!< List of population analyses
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficient matrix
    const ab_matrix &m_dm0; //!< Ground state density matrix
    const ab_matrix &m_dm; //!< Excited state (or difference) density matrix
    const ev_printer_i &m_prno; //!< Formating object of NO summary
    const ev_printer_i &m_prndo; //!< Formating object of NDO summary
    bool m_is_diff; //!< Density matrix is difference density

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param dm0 Ground state density matrix
        \param dm State or difference density matrix
        \param prno NO summary printer
        \param prndo NDO summary printer
        \param is_diff Is difference density?
     **/
    analyse_esdm(
        const arma::Mat<double> &s,
        const ab_matrix &c,
        const ab_matrix &dm0,
        const ab_matrix &dm,
        const ev_printer_i &prno,
        const ev_printer_i &prndo,
        bool is_diff = true);

    /** \brief Register population analyses that should be performed
        \param name Name for population analysis
        \param pa Population analysis
        \param pr Population printer / formatter
        \param fl Flag for which density the population analysis is used
     **/
    void do_register(const std::string &name, const pop_analysis_i &pa,
        const pop_printer_i &pr, pa_flag fl);

    /** \brief Performs density matrix analyses
        \param pr Density and NTO export / printer
        \param out Output stream

        Perform the following analyses:
        - NO analysis
        - NDO analysis
        - Population analysis
     **/
    void perform(export_data_i &pr, std::ostream &out) const;
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_SDM_H
