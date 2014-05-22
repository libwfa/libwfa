#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

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
    ctnumbers m_ct; //!< CT number analysis
    nto_analysis m_nto; //!< NTO analysis
    export_densities_i &m_pr_d; //!< Density printer
    export_orbitals_i &m_pr_o; //!< Orbital printer

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
        const ctnum_analysis_i &ctnum, export_densities_i &pr_d,
        export_orbitals_i &pr_o, ev_printer_i &pr_nto,
        ctnum_printer_i &pr_ct);

    /** \brief Performs transition density matrix analyses
        \param tdm Transition density matrix
        \param[out] av Average particle and hole density matrices (particle first)

        Perform the following analyses:
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
        - Export of TDM, EDM, and HDM
        - EDM and HDM are added to av
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &av);

    /** \brief Perform CT number analysis
        \param tdm Transition density matrix
        \param pr Printer of CT number data
     **/
    void analyse_ctnum(const ab_matrix &tdm, std::ostream &out) const {
        m_ct.perform(tdm, out);
    }

    /** \brief Performs NTO analysis
        \param tdm Transition density matrix
        \param pr_d Density matrix export / print
        \param pr_o NTO export / print
        \param out Output stream

        EDM and HDM are exported and discarded afterwards.
     **/
    void analyse_nto(const ab_matrix &tdm, export_densities_i &pr_d,
        export_orbitals_i &pr_o, std::ostream &out) const;

};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
