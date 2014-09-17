#ifndef SANTO_ANALYSIS_H
#define SANTO_ANALYSIS_H

#include "santo_analysis.h"
#include <libwfa/export/ev_printer_i.h>
#include <libwfa/export/export_data_i.h>

namespace libwfa {

/** \brief Class for transformation of a TDM with respect to SA-NTOs
 *
 * Transformation of the transition density matrix with respect to
 *  state-averaged natural transition orbitals.
 *
 * \ingroup libwfa
 */
class santo_analysis {
private:
    ab_matrix m_uinv, m_vinv_t; //!< Hole and electron transformation matrices

public:
    /** \brief constructor
     * 
     *  Read in the data, diagonalize the state-averaged hole and electron
     *   density matrices, and construct the transformation matrices.
     *
     *  If edm and hdm contain different alpha and beta components, then
     *   separate alpha and beta SA-NTOs are created.
     *  If, on the other hand, spin-averaged SA-NTOs are wanted,
     *   spin-averaged edm and hdm matrices have to be passed.
     * 
     *  \param s Overlap matrix
     *  \param c Orbital coefficient matrix for transform in orthogonal basis
     *  \param edm State-averaged electron density matrix
     *  \param hdm State-averaged hole density matrix
     *  \param evpr Formatting object for eigenvalue print-out
     *  \param dpr Printer of SA-NTOs
     *  \param out Output stream
     **/
    santo_analysis(const arma::mat &s, const ab_matrix &c,
         const ab_matrix &edm, const ab_matrix &hdm,
         const ev_printer_i &evpr, export_data_i &opr,
         std::ostream &out);
    
    void stave_nto_header(std::ostream& out);
    
    /** \brief Analyze a transition density matrix
     * 
     *  Call decompose and print routines
     * 
     *  \param tdm Transition density matrix
     *  \param out Output stream
     */
    void perform(const ab_matrix &tdm, std::ostream &out) {
        ab_matrix x(tdm.is_alpha_eq_beta());
        decompose(tdm, x);
        print(x, out);
    }
    
    /** \brief Return the SA-NTO decomposition
     *
     * \param tdm[in] transition density matrix to be decomposed
     * \param x[out]  decomposed tdm
     *
     */
    void decompose(const ab_matrix &tdm, ab_matrix &x) {
        x = m_uinv * tdm * m_vinv_t;
    }
    
    /** \brief Return a transformation matrix
     * 
     *  \param elec electron (true) or hole (false) matrix
     **/
    ab_matrix &get_trans(bool elec) {
        if (elec) return m_vinv_t;
        else return m_uinv;
    }
    
    /** \brief Print results of SA-NTO decomposition
     *
     *  \param x   decomposed tdm
     *  \param out Output stream
     */
    void print(const ab_matrix &x, std::ostream &out);
    
private:
    /** \brief Print results of SA-NTO decomposition (alpha or beta)
     *
     *  \param xm   decomposed tdm (alpha or beta part)
     *  \param out  Output stream
     *  \param dfac Set to 2. in the case of spin-traced NTOs
     */
    void print(const arma::mat &xm, std::ostream &out, double dfac = 1.);
    
};

} // namespace libwfa


#endif /* SANTO_ANALYSIS_H_ */
