#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

#include "no_analysis.h"
#include "ndo_analysis.h"
#include "pop_analysis_i.h"
#include "pop_printer_i.h"

namespace libwfa {

/** \brief Combines various analyses of a state density matrix

    \ingroup libwfa
 **/
class analyse_sdm {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;

private:
    const pop_analysis_i &m_pop; //!< Population analysis to perform
    ndo_analysis m_ndo; //!< NDO analysis
    no_analysis m_no; //!< NO analysis
    const ab_matrix &m_gs_dm; //!< Ground state density matrix

public:
    analyse_sdm(const pop_analysis_i &pop, const arma::Mat<double> &s,
        const ab_matrix &c, const ab_matrix &gs_dm) :
        m_pop(pop), m_no(c), m_ndo(s, c), m_gs_dm(gs_dm) { }

    /** \brief Performs density matrix analyses
        \param sdm State density matrix
        \param pr_d Density matrix export / print
        \param pr_o Orbital export / print
        \param pr_e1 Printer of NO summary
        \param pr_e2 Printer of NDO summary
        \param pr_p Printer of population data

        Perform the following analyses:
        - Export of TDM
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
     **/
    void perform(const ab_matrix &sdm,
        export_densities_i &pr_d, export_orbitals_i &pr_o,
        ev_printer_i &pr_e1, ev_printer_i &pr_e2, pop_printer_i &pr_p) const;

    /** \brief Perform population analysis
        \param sdm State density matrix
        \param ad Attachment / detachment densities
        \param pr Printer of population data
     **/
    void pop_analysis(const ab_matrix &sdm, const ab_matrix_pair &ad,
            pop_printer_i &pr) const;

    /** \brief Performs NO analysis
        \param sdm State density matrix
        \param pr_o NTO export / print
        \param pr_e Printer of NTO summary
     **/
    void no_analysis(const ab_matrix &sdm, export_orbitals_i &pr_o,
        ev_printer_i &pr_e) const {

        m_no.perform(sdm, pr_o, pr_e);
    }

    /** \brief Performs NO analysis
        \param sdm State density matrix
        \param ad Attachment / detachment density matrices
        \param ndo_print NDO export / print
        \param pr Printer of NDO summary
     **/
    void ndo_analysis(const ab_matrix &sdm, ab_matrix_pair &ad,
            export_orbitals_i &pr_o, ev_printer_i &pr_e) const {

        ab_matrix ddm(sdm);
        ddm -= m_gs_dm;
        m_ndo.perform(ddm, ad, pr_o, pr_e);
    }
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
