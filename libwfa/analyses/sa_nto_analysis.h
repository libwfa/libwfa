#ifndef LIBWFA_SA_NTO_ANALYSIS_H
#define LIBWFA_SA_NTO_ANALYSIS_H

#include "nto_analysis.h"

namespace libwfa {

/** \brief Class for transformation of a TDM with respect to SA-NTOs

    Transformation of the transition density matrix with respect to
    state-averaged natural transition orbitals.

    \ingroup libwfa
 **/
class sa_nto_analysis {
private:
    ab_matrix m_ul, m_ur; //!< Hole and electron transformation matrices

public:
    /** \brief Constructor
        \param s Overlap
        \param nto NTO data

        Construct the NTO data to form transformation matrices.
     **/
    sa_nto_analysis(const arma::mat &s, const nto_analysis &ntos);

    /** \brief Analyze a transition density matrix
        \param tdm Transition density matrix
        \param out Output stream
        \param thresh Threshold

        Decomposes the transition density matrix and analyse the result
     **/
    void analyse(const ab_matrix &tdm, std::ostream &out,
        double thresh = 1e-2) const;
    
private:
    /** \brief Analyse decomposed transition density matrices)
        \param out Output stream
        \param x Decomposed TDM (alpha or beta part)
        \param c Scaling factor
        \param thresh Threshold of important contributions
     **/
    static void analysis(std::ostream &out, const arma::mat &x,
        double c, double thresh);
    
};

} // namespace libwfa


#endif /* SANTO_ANALYSIS_H_ */
