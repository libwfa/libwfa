#ifndef LIBWFA_ANALYSE_DDM_H
#define LIBWFA_ANALYSE_DDM_H

#include <map>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/export/pop_printer_i.h>

namespace libwfa {

/** \brief Combines various analyses of a difference density matrix

    \ingroup libwfa
 **/
class analyse_ddm {
private:
    struct pa {
        const pop_analysis_i &analysis;
        const pop_printer_i &printer;

        pa(const pop_analysis_i &a, const pop_printer_i &p) :
            analysis(a), printer(p) { }
    };
    typedef std::map<std::string, pa> pa_map_t;

private:
    pa_map_t m_lst; //!< List of population analyses
    const arma::Mat<double> &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficient matrix
    const ab_matrix &m_ddm; //!< Difference density matrix
    const ev_printer_i &m_prndo; //!< Formating object of NDO summary

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c MO coefficients
        \param ddm Difference density matrix
        \param pr NDO summary formating object
     **/
    analyse_ddm(
        const arma::Mat<double> &s,
        const ab_matrix &c,
        const ab_matrix &ddm,
        const ev_printer_i &pr);

    /** \brief Register population analyses to be performed
        \param name Name for population analysis
        \param pa Population analysis
        \param pr Population printer / formatter
     **/
    void do_register(const std::string &name,
        const pop_analysis_i &pa, const pop_printer_i &pr);

    /** \brief Performs density matrix analyses
        \param pr Density and orbital export / printer
        \param out Output stream

        Perform the following analyses:
        - NDO analysis
        - Population analysis
     **/
    void perform(export_data_i &pr, std::ostream &out) const;
};

} // namespace libwfa

#endif // LIBWFA_ANALYSE_SDM_H
