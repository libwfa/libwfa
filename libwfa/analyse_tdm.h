#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

#include <map>
#include <libwfa/analyses/ctnumbers.h>
#include <libwfa/analyses/nto_analysis.h>
#include <libwfa/export/export_densities_i.h>

namespace libwfa {

/** \brief Combines various transition density matrix analyses

    \ingroup libwfa
 **/
class analyse_tdm {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    struct cna {
        const ctnum_analysis_i &analysis;
        const ctnum_printer_i &printer;

        cna(const ctnum_analysis_i &a, const ctnum_printer_i &p) :
            analysis(a), printer(p) { }
    };
    typedef std::map<std::string, cna> cna_map_t;


private:
    cna_map_t m_lst; //!< List of CT number analysis
    const arma::Mat<double> &m_s;
    nto_analysis m_nto; //!< NTO analysis

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Coefficient matrix
        \param ctnum Charge transfer number analysis
        \param pr_d Density printer
        \param pr_o Orbital printer
        \param pr_nto NTO summary printer
        \param pr_ct CT number printer
     **/
    analyse_tdm(const arma::Mat<double> &s, const ab_matrix &c,
        ev_printer_i &pr_nto);

    /** \brief Register CT number analysis to be performed
        \param name Name for CT number analysis
        \param ana CT number analysis
        \param pr CT number printer
     **/
    void do_register(const std::string &name, const ctnum_analysis_i &ana,
        const ctnum_printer_i &pr);

    /** \brief Performs transition density matrix analyses
        \param tdm Transition density matrix
        \param[out] av Average particle and hole density matrices (particle first)

        Perform the following analyses:
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
        - Export of TDM, EDM, and HDM
        - EDM and HDM are added to av
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &dpr, export_orbitals_i &opr, std::ostream &out);

    /** \brief Perform CT number analysis
        \param tdm Transition density matrix
        \param pr Printer of CT number data
     **/
    void analyse_ctnum(const ab_matrix &tdm, std::ostream &out) const;

    /** \brief Performs NTO analysis
        \param tdm Transition density matrix
        \param dpr Density matrix export / print
        \param opr NTO export / print
        \param out Output stream

        EDM and HDM are exported and discarded afterwards.
     **/
    void analyse_nto(const ab_matrix &tdm, export_densities_i &dpr,
        export_orbitals_i &opr, std::ostream &out) const;

};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
