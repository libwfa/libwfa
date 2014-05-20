#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

#include "ctnumbers.h"
#include "export_densities_i.h"
#include "nto_analysis.h"

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

public:
    analyse_tdm(const ctnum_analysis_i &ctnum,
        const arma::Mat<double> &s, const ab_matrix &c) :
        m_ct(ctnum, s), m_nto(s, c) { }

    /** \brief Performs transition density matrix analyses
        \param tdm Transition density matrix
        \param av Average particle and hole density matrices (particle first)
        \param pr_d Density matrix export / print
        \param pr_o NTO export / print
        \param pr_e Printer of NTO summary
        \param pr_c Printer of CT number data

        Perform the following analyses:
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
        - Export of TDM, EDM, and HDM
        - EDM and HDM are added to av
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &pr_d, export_orbitals_i &pr_o,
        ev_printer_i &pr_e, ctnum_printer_i &pr_c) const;

    /** \brief Perform CT number analysis
        \param tdm Transition density matrix
        \param pr Printer of CT number data
     **/
    void ctnumbers(const ab_matrix &tdm, ctnum_printer_i &pr) const {
        m_ct.perform(tdm, pr);
    }

    /** \brief Performs NTO analysis
        \param tdm Transition density matrix
        \param pr_d Density matrix export / print
        \param pr_o NTO export / print
        \param pr_e Printer of NTO summary

        EDM and HDM are exported and discarded afterwards.
     **/
    void nto_analysis(const ab_matrix &tdm, export_densities_i &pr_d,
        export_orbitals_i &pr_o, ev_printer_i &pr_e) const;

};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
