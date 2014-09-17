#ifndef LIBWFA_NO_ANALYSIS_H
#define LIBWFA_NO_ANALYSIS_H

#include <libwfa/core/transformations_dm.h>
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>

namespace libwfa {


/** \brief Perform complete NO analysis of a state density matrix

     \ingroup libwfa
 **/
class no_analysis {
private:
/*    const arma::mat &m_s; //!< Overlap matrix
    const ab_matrix &m_c; //!< MO coefficients
    const ab_matrix &m_sdm; //!< State density matrix
    const ev_printer_i &m_pr; //!< Formating object*/
    ab_vector m_e; //!< Eigenvalues, i.e. occupation numbers
    ab_matrix m_u; //!< Eigenvectors, i.e. natural orbitals
    ab_vector m_est; //!< spin-traced eigenvalues
    ab_matrix m_ust; //!< spin-traced eigenvectors
    


public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Orbital coefficient matrix for transform in orthogonal basis
        \param sdm State density matrix
        \param pr Formatting object
     **/
    no_analysis(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &sdm) :
        m_e(c.is_alpha_eq_beta()), m_u(c.is_alpha_eq_beta()),
        m_est(c.is_alpha_eq_beta()), m_ust(c.is_alpha_eq_beta()) {
            diagonalize_dm(s, c, sdm, m_e, m_u);
            
            // Perform spin-traced calculation
            if (! c.is_alpha_eq_beta()) {
                ab_matrix sdm2(true);
                sdm2.alpha() = sdm.alpha() + sdm.beta();

                diagonalize_dm(s, c, sdm2, m_est, m_ust);

                m_est.alpha() *= 0.5;
            }
            
        }

    /** \brief Return eigenvalues
     * 
     *  \param spint return of spin-traced results
     **/
    ab_vector &get_eigval(bool spint) {
        if (spint && ! m_e.is_alpha_eq_beta()) return m_est;
        else return m_e;
    }
        
    /** \brief Return eigenvectors
     * 
     *  \param spint return of spin-traced results
     **/
    ab_matrix &get_eigvect(bool spint) {
        if (spint && ! m_e.is_alpha_eq_beta()) return m_ust;
        else return m_u;
    }
    
    /** \brief Perform NTO analysis
        \param evpr Eigenvalue printer
        \param opr Printer of NOs
        \param out Output stream
     **/
    void perform(const ev_printer_i &evpr, export_data_i &opr, 
                 std::ostream &out) const;
};


} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_H
