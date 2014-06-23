#ifndef LIBWFA_ANALYSE_OPTDM_H
#define LIBWFA_ANALYSE_OPTDM_H

#include <map>
#include <libwfa/analyses/ctnum_analysis_i.h>
#include <libwfa/export/ctnum_printer_i.h>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>
//ex_analyse part
#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/export/ex_ana_printer.h>

namespace libwfa {

/** \brief Combines various transition density matrix analysis

    \ingroup libwfa
 **/
class analyse_optdm {
private:
    struct cna {
        const ctnum_analysis_i &analysis;
        const ctnum_printer_i &printer;

        cna(const ctnum_analysis_i &a, const ctnum_printer_i &p) :
            analysis(a), printer(p) { }
    };
    typedef std::map<std::string, cna> cna_map_t;


private:
    cna_map_t m_ca; //!< List of CT number analysis
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficient matrix
    const ab_matrix &m_tdm; //!< Transition density matrix
    const ev_printer_i *m_nto; //!< Formating object of NTO summary


public:
    /** \brief Constructor for additional exciton analysis
        \param s Overlap Matrix
        \param c Coefficient matrix
        \param tdm Transistion density matrix
     */
    analyse_optdm(const arma::Mat<double> &s, const ab_matrix &c,
            const ab_matrix &tdm);

    /** \brief Register NTO analysis
        \param pr NTO summary printer
     **/
    void do_register(const ev_printer_i &pr) { m_nto = &pr; }

    /** \brief Register CT number analysis to be performed
        \param name Name for CT number analysis
        \param ana CT number analysis
        \param pr CT number printer
     **/
    void do_register(const std::string &name, const ctnum_analysis_i &ana,
        const ctnum_printer_i &pr);
    /** \brief Performs transition density matrix analyses
        \param edm_av Average electron density matrix
        \param hdm_av Average hole density matrix
        \param pr Export / printer for densities and orbitals
        \param con Contract Interface
        \param out Output stream

        Perform the following analyses:
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
        - Export of TDM, EDM, and HDM
        - EDM and HDM are added to average density matrices
        - Exciton analysis
     **/
    void perform(ab_matrix &edm_av, ab_matrix &hdm_av,
        export_data_i &pr, const multipol_con_i &con, std::ostream &out);

    /** \brief Performs transition density matrix analyses
        \param pr Export / printer for densities and orbitals
        \param con Contract Interface
        \param out Output stream

        Perform the following analyses:
        - NTO analysis (\sa nto_analysis.h
        - CT number analysis (\sa ctnumbers.h)
        - Export of TDM, EDM, and HDM
        - Exciton analysis
     **/
    void perform(export_data_i &pr, const multipol_con_i &con,
            std::ostream &out);
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_OPTDM_H
