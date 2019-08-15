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


#ifndef LIBWFA_WF_ANALYSIS_H
#define LIBWFA_WF_ANALYSIS_H

#include <unordered_map>
#include <map>
#include <libwfa/analyses/sa_nto_analysis.h>
#include "wf_analysis_data_i.h"

namespace libwfa {


/** \brief Wave function analysis class

    Wrapper class to perform any or all available analyses. The types of
    analyses that are performed is determined by the analysis data object.

    \ingroup libwfa
 **/
class wf_analysis {
private:
    typedef wf_analysis_data_i::orbital_params orbital_params;

private:
    std::unique_ptr<wf_analysis_data_i> m_h; //!< Analysis data
    std::unique_ptr<sa_nto_analysis> m_sa; //!< State-averaged NTO analysis
    ab_matrix m_edm_av; //!< Averaged electron density
    ab_matrix m_hdm_av; //!< Averaged hole density
    bool m_init_av; //!< Whether the above are initialized?

    /** \brief Data relating to OmFrag analysis
     **/
    struct frag_data {
        std::string state_name;
        double dE_eV, f = 0.0;
        double om_tot;
        arma::mat om_frag;
        std::unordered_map<std::string, double> descriptor;
    };
    std::unordered_map<int, std::map<std::string, frag_data> > frag_data_all; //!< final output to be printed


public:
    /** \brief Initializes the wave function analysis
        \param h Analysis data object
     **/
    wf_analysis(wf_analysis_data_i *h) : m_h(h), m_sa(), m_init_av(false) { }

    /** \brief Perform analysis of state and difference density matrix
        \param out Output stream
        \param name Name of state (useable as filename)
        \param desc Description of state (one-line comment)
        \param ddm Difference density matrix in AO
        \param dm0 Ground state density matrix in AO
     **/
    void analyse_opdm(std::ostream &out,
        const std::string &name, const std::string &desc,
        const ab_matrix &ddm, const ab_matrix &dm0);

    /** \brief Perform analysis of state density matrix
            \param out Output stream
            \param name Name of state (useable as filename)
            \param desc Description of state (one-line comment)
            \param sdm State density matrix in AO
     **/
    void analyse_opdm(std::ostream &out, const std::string &name,
        const std::string &desc, const ab_matrix &sdm);

    /** \brief Perform analysis of transition density matrix
        \param out Output stream
        \param name Name of state (useable as filename)
        \param desc Description of state (one-line comment)
        \param tdm Transition density matrix in AO
        \param energy Excitation energy (eV)
     **/
    void analyse_optdm(std::ostream &out, const std::string &name,
        const std::string &desc, const ab_matrix &tdm, double energy=0.);

    /** \brief Constructs state-averaged NTOs and sets up the analysis
        \param out Output stream
        \return True, if successful
     **/
    bool setup_sa_ntos(std::ostream &out);

    /** \brief Perform SA-NTO analysis of transition density matrix
        \param out Output stream
        \param tdm Transition density matrix in AO
        \return True, if successful
     **/
    bool post_process_optdm(std::ostream &out, const ab_matrix &tdm);

    /** \brief Reset to original state (at construction)
     **/
    void reset() {
        m_init_av = false;
        delete m_sa.release();
    }

    /** \brief Print a TheoDORE-style summary of the CT number analysis
     **/
     void print_summary(std::ostream &out, const int &prec=6, const int &width=10);

     /** \brief Print an OmFrag.txt file as understood by TheoDORE
      **/
      void print_om_frag(std::ostream &out, const std::string ofile="OmFrag.txt");

private:
    void add_to_average(const ab_matrix &edm, const ab_matrix &hdm);
};


/** \brief Holder for static analysis object

    The analysis object needs to be initialized, for the C-style interface to
    work properly.
 **/
struct wf_analysis_static {
    static std::unique_ptr<wf_analysis> analysis; //!< Static analysis object
};


} // namespace libwfa

#endif // LIBWFA_WF_ANALYSIS_H
