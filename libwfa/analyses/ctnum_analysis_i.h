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


#ifndef LIBWFA_CTNUM_ANALYSIS_I_H
#define LIBWFA_CTNUM_ANALYSIS_I_H

#include <armadillo>
#include <cstdlib>
#include <unordered_map>

namespace libwfa {

/** \brief Base class for charge transfer (CT) number analysis.

    \ingroup libwfa
 **/
class ctnum_analysis_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~ctnum_analysis_i() { }
    
    /** \brief Return the dimension of the CT number matrix (e.g. # atoms)
     **/
    virtual size_t size() const = 0;

    /** \brief Perform CT number analysis for given transition density matrix
        \param[in] tdm Transition density matrix in AO basis
        \param[out] om Result CT number matrix

        Perform the CT number analysis for the supplied transition density
        matrix in AO basis. Any implementation should adjust the output matrix
        to the correct size.
     **/
    virtual void perform(const arma::mat &tdm, arma::mat &om) const = 0;

    virtual std::unordered_map<std::string, double> compute_desc(const std::vector<double> &om_tot,
                                                                const arma::mat &om) const = 0;
};

} // namespace libwfa

#endif
