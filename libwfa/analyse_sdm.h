#ifndef LIBWFA_ANALYSE_SDM_H
#define LIBWFA_ANALYSE_SDM_H

#include <map>
#include <libwfa/analyses/no_analysis.h>
#include <libwfa/analyses/ndo_analysis.h>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/export/pop_printer_i.h>

namespace libwfa {

/** \brief Combines various analyses of state density matrices

    \ingroup libwfa
 **/
class analyse_sdm {
public:
    typedef std::pair<ab_matrix, ab_matrix> ab_matrix_pair;
    typedef enum {
        pa_dm = 1, pa_ad = 2, both = 3
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
    pa_map_t m_lst; //!< List of population analysis
    const ab_matrix &m_dm0; //!< Ground state density matrix
    no_analysis m_no; //!< NO analysis
    ndo_analysis m_ndo; //!< NDO analysis

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param dm0 Ground state density matrix
        \param pr_no NO summary printer
        \param pr_ndo NDO summary printer
     **/
    analyse_sdm(const arma::Mat<double> &s, const ab_matrix &c,
        const ab_matrix &dm0, const ev_printer_i &pr_no,
        const ev_printer_i &pr_ndo);

    /** \brief Register population analysis to be performed
        \param name Name for population analysis
        \param pa Population analysis
        \param pr Population printer
        \param fl Flag to apply
     **/
    void do_register(const std::string &name, const pop_analysis_i &pa,
            const pop_printer_i &pr, pa_flag fl);

    /** \brief Performs density matrix analyses
        \param[in] dm State or difference density matrix
        \param[out] dpr Density export / printer
        \param[out] opr NTO export / printer
        \param[out] is_diff Is the first argument a difference density matrix?

        Perform the following analyses:
        - NO analysis
        - NDO analysis
        - Population analysis
     **/
    void perform(const ab_matrix &dm, export_densities_i &dpr,
            export_orbitals_i &opr, std::ostream &out, bool is_diff = true);

    /** \brief Perform population analysis
        \param[in] sdm State density matrix
        \param[out] ad Attachment / detachment densities
        \param[out] out Output stream
     **/
    void analyse_pop(const ab_matrix &sdm, const ab_matrix_pair &ad,
            std::ostream &out) const;

    /** \brief Performs NO analysis
        \param[in] sdm State density matrix
        \param[out] opr NTO export / print
        \param[out] out Output stream
     **/
    void analyse_no(const ab_matrix &sdm, export_orbitals_i &opr,
        std::ostream &out) const {

        m_no.perform(sdm, opr, out);
    }

    /** \brief Performs NDO analysis
        \param[in] ddm Difference density matrix
        \param[out] ad Attachment / detachment density matrices
        \param[out] opr NDO export / print
        \param[out] out Output stream
     **/
    void analyse_ndo(const ab_matrix &ddm, ab_matrix_pair &ad,
        export_orbitals_i &opr, std::ostream &out) const {

        m_ndo.perform(ddm, ad, opr, out);
    }
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_SDM_H
