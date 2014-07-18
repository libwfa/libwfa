#ifndef SANTO_ANALYSIS_H
#define SANTO_ANALYSIS_H

#include "nto_analysis.h"
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
    ab_matrix m_uinv_t, m_vinv; //!< Hole and electron transformation matrices

public:
    /** \brief constructor
     * 
     *  Read in the data, diagonalize the state-averaged hole and electron
     *   density matrices, and construct the transformation matrices.
     * 
     *  \param s Overlap matrix
     *  \param c Orbital coefficient matrix for transform in orthogonal basis
     *  \param edm State-averaged electron density matrix
     *  \param edm State-averaged hole density matrix
     *  \param evpr Formatting object for eigenvalue print-out
     *  \param dpr Printer of SA-NTOs
     *  \param out Output stream
     **/
    santo_analysis(const arma::Mat<double> &s, const ab_matrix &c,
         const ab_matrix &edm, const ab_matrix &hdm,
         const ev_printer_i &evpr, export_data_i &opr,
         std::ostream &out);
    
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
     * \param tdm transition density matrix to be decomposed
     * \param x   decomposed tdm
     *
     */
    void decompose(const ab_matrix &tdm, ab_matrix &x) {
        x = m_uinv_t * tdm * m_vinv;
    }
    
    /** \brief Print results of SA-NTO decomposition
     *
     *  \param x   decomposed tdm
     *  \param out Output stream
     */
    void print(const ab_matrix &x, std::ostream &out);
};

} // namespace libwfa


#endif /* SANTO_ANALYSIS_H_ */
