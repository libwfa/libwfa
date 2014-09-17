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
        \param nto State-averaged NTO data

        Construct the NTO data to form transformation matrices.
     **/
    sa_nto_analysis(const arma::mat &s, const nto_analysis &ntos);

    /** \brief Retrieve left transformation matrix
     **/
    const ab_matrix &get_transf_l() { return m_ul; }

    /** \brief Retrieve right transformation matrix
     **/
    const ab_matrix &get_transf_r() { return m_ur; }

    /** \brief Decomposes a transition density matrix
        \param tdm Transition density matrix
        \param xdm Decomposed matrix
     **/
    void decompose(const ab_matrix &tdm, ab_matrix &xdm) const;

    /** \brief Analyse a transition density matrix
        \param out Output stream
        \param tdm Transition density matrix
        \param xdm Decomposed matrix
        \param thresh Threshold

        Decomposes the transition density matrix and analyses the result
     **/
    void analyse(std::ostream &out, const ab_matrix &tdm,
        ab_matrix &xdm, double thresh = 1e-2) const;

    /** \brief Analyse a transition density matrix
        \param out Output stream
        \param tdm Transition density matrix
        \param thresh Threshold

        Decomposes the transition density matrix and analyses the result
     **/
    void analyse(std::ostream &out, const ab_matrix &tdm,
        double thresh = 1e-2) const {

        ab_matrix xdm;
        analyse(out, tdm, xdm, thresh);
    }

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


#endif // LIBWFA_SA_NTO_ANALYSIS_H
