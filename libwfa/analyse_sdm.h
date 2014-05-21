#ifndef LIBWFA_ANALYSE_SDM_H
#define LIBWFA_ANALYSE_SDM_H

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
    const ab_matrix &m_dm0; //!< Ground state density matrix
    ndo_analysis m_ndo; //!< NDO analysis
    no_analysis m_no; //!< NO analysis
    const pop_analysis_i &m_pop; //!< Population analysis to perform
    export_densities_i &m_pr_d; //!< Density printer
    export_orbitals_i &m_pr_o; //!< Orbital printer
    ev_printer_i &m_pr_no; //!< NO summary printer
    ev_printer_i &m_pr_ndo; //!< NDO summary printer
    pop_printer_i &m_pr_p; //!< Printer of population data

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param dm0 Ground state density matrix
        \param pop Population analysis class
        \param pr_d Density printer
        \param pr_o Orbital printer
        \param pr_no NO summary printer
        \param pr_nto NTO summary printer
        \param pr_p Printer of population data
     **/
    analyse_sdm(
        const arma::Mat<double> &s,
        const ab_matrix &c,
        const ab_matrix &dm0,
        const pop_analysis_i &pop,
        export_densities_i &pr_d,
        export_orbitals_i &pr_o,
        ev_printer_i &pr_no,
        ev_printer_i &pr_nto,
        pop_printer_i &pr_p);

    /** \brief Performs density matrix analyses
        \param dm State or difference density matrix
        \param is_diff Is difference density matrix?

        Perform the following analyses:
        - NO analysis
        - NDO analysis
        - Population analysis
     **/
    void perform(const ab_matrix &dm, bool is_diff = true);

    /** \brief Perform population analysis
        \param sdm State density matrix
        \param ad Attachment / detachment densities
        \param pr Printer of population data
     **/
    void analyse_pop(const ab_matrix &sdm,
        const ab_matrix_pair &ad, pop_printer_i &pr) const;

    /** \brief Performs NO analysis
        \param sdm State density matrix
        \param pr_o NTO export / print
        \param pr_e Printer of NTO summary
     **/
    void analyse_no(const ab_matrix &sdm, export_orbitals_i &pr_o,
        ev_printer_i &pr_e) const {

        m_no.perform(sdm, pr_o, pr_e);
    }

    /** \brief Performs NDO analysis
        \param ddm Difference density matrix
        \param[out] ad Attachment / detachment density matrices
        \param pr_o NDO export / print
        \param pr_e Printer of NDO summary
     **/
    void analyse_ndo(const ab_matrix &ddm, ab_matrix_pair &ad,
        export_orbitals_i &pr_o, ev_printer_i &pr_e) const {

        m_ndo.perform(ddm, ad, pr_o, pr_e);
    }
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_SDM_H
