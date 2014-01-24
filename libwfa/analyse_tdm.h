#ifndef LIBWFA_ANALYSE_TDM_H
#define LIBWFA_ANALYSE_TDM_H

#include "ctnumbers.h"
#include "nto_analysis.h"

namespace libwfa {

/** \brief Combines various transition density matrix analyses

    Performs the analyses of the current transition density according to the
    parameters.

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
        \param av Average electron and hole density matrices
        \param dm_print Density matrix export / print
        \param nto_print NTO export / print
        \param prn Printer of NTO summary
        \param prct Printer of CT number data

        Perform the following analyses:
        - Export of TDM
        - NTO analysis (\sa nto_analysis.h)
        - CT number analysis (\sa ctnumbers.h)
     **/
    void perform(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &dm_print, export_orbitals_i &nto_print,
        nto_data_i &prn, ctnum_data_i &prct) const;

    /** \brief Perform CT number analysis
        \param tdm Transition density matrix
        \param pr Printer of CT number data
     **/
    void ctnumbers(const ab_matrix &tdm, ctnum_data_i &pr) const {
        m_ct.perform(tdm, pr);
    }

    /** \brief Performs NTO analysis
        \param tdm Transition density matrix
        \param av Average electron and hole density matrices
        \param dm_print Density matrix export / print
        \param nto_print NTO export / print
        \param prn Printer of NTO summary
     **/
    void nto_analysis(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &dm_print, export_orbitals_i &nto_print,
        nto_data_i &pr) const {
        m_nto.perform(tdm, av, dm_print, nto_print, pr);
    }
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_TDM_H
