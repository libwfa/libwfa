//************************************************************************
//* This file is part of libwfa.                                         *
//*                                                                      *
//* libwfa is free software; you can redistribute and/or modify          *
//* it under the terms of the BSD 3-Clause license.                      *
//* libwfa is distributed in the hope that it will be useful, but it     *
//* is provided "as is" and without any express or implied warranties.   *
//* For more details see the full text of the license in the file        *
//* LICENSE.                                                             *
//*                                                                      *
//* Copyright (c) 2014, F. Plasser and M. Wormit. All rights reserved.   *
//* Modifications copyright (C) 2019, Loughborough University.           *
//************************************************************************


#ifndef LIBWFA_CTNUM_ANALYSIS_H
#define LIBWFA_CTNUM_ANALYSIS_H

#include <vector>
#include <unordered_map>
#include "ctnum_analysis_i.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Basic implementation of CT number analysis. Provided a mapping of basis
    functions to molecular fragments (e.g. atoms) the class computes for every
    \f$ \Omega \f$ matrix in AOs the CT number matrix in terms of the fragments.

    \ingroup libwfa
 **/
class ctnum_analysis : public ctnum_analysis_i {
private:
    typedef std::vector<std::vector<int>> ivector;

    size_t m_nparts; //!< Number of parts
    const arma::mat &m_s; //!< Overlap matrix
    arma::mat m_s_sqrt; //!< square root of overlap matrix
    const arma::uvec &m_b2p; //!< Map of basis functions to fragments
    const std::string ctnum_method; //!< Formula used to calculate CT number

    const std::vector<std::string> prop_list;
    const ivector at_lists;


public:
    /** \brief Constructor
        \param s Overlap matrix
        \param b2p Map of basis functions to molecular parts or fragments
        \param method Formula used to calculate CT number
        \param prop_list Properties to be computed
        \param at_lists Fragment definitions
     **/
    ctnum_analysis(const arma::mat &s, const arma::uvec &b2p, const std::string &method,
            const std::vector<std::string> &prop_list, const std::vector<std::vector<int>> &at_lists);

    /** \brief Constructor
        \param s Overlap matrix
        \param b2p Map of basis functions to molecular parts or fragments
        \param method Formula used to calculate CT number
     **/
    ctnum_analysis(const arma::mat &s, const arma::uvec &b2p, const std::string &method);

    /** \brief Virtual destructor
     **/
    virtual ~ctnum_analysis() { }

    /** \copydoc ctnum_analysis_i::size
     **/
    virtual size_t size() const { return m_nparts; }

    /** \copydoc ctnum_analysis_i::perform
     **/
    virtual void perform(const arma::mat& tdm, arma::mat &om, double &Phe) const;

    /** \brief Summation of the Omega matrix over fragments
        \param[in] om_at Omega matrix over atoms
     **/
    virtual arma::mat compute_omFrag(const arma::mat &om_at) const;

    /** \brief Compute general Frobenius scalar product
        \param[in] D1 input matrix
        \param[in] D2 input matrix

        Output: tr(D1 * S * D2.T * S)
    **/
    virtual double compute_DSDS(const arma::mat &D1, const arma::mat &D2) const;

    /** \brief Number of fragments
     **/
    virtual size_t n_frags() const { return at_lists.size(); }

    /** \brief Compute the TheoDORE-style descriptors
     **/
   virtual std::unordered_map<std::string, double> compute_descriptors(const double om_tot,
                                                               const arma::mat &om_frag) const;

//private:
    /** \brief Forms omega matrix from a transition density matrix
        \param[in] s Overlap matrix
        \param[in] s_sqrt Square root of overlap matrix
        \param[in] tdm Transition density matrix
        \param[in] method Formula used to calculate CT number
        \param[out] om Omega matrix

        The function implements the transforms of the transition density
        matrix into the \f$\Omega\f$ matrix, which is further used to compute
        the charge transfer numbers.

        The implementation uses the new formula
        \f[
        \Omega_{\mu\nu} = 0.5 \left[
            (\mathbf{D}\mathbf{S})_{\mu\nu} \times (\mathbf{S}\mathbf{D})_{\mu\nu}
            + D_{\mu\nu} \times (\mathbf{S}\mathbf{D}\mathbf{S})_{\mu\nu}
            \right]
        \f]
         from [JCP(2014), DOI: 10.1063/1.4885819] rather than the original
         formula from [JCTC(2012), 8, 2777].

         In addition a Lowdin orthogonalization is implemented.

        The output matrices are reshaped and resized as required.
     **/
    static void form_om(const arma::mat &s, const arma::mat &s_sqrt,
            const arma::mat &tdm, const std::string &method, arma::mat &om);
};

} // namespace libwfa

#endif
