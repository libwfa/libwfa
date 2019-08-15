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
//* Copyright (C) 2019, Loughborough University                          *
//************************************************************************

#ifndef LIBWFA_OM_DESCRIPTOR_H
#define LIBWFA_OM_DESCRIPTOR_H

#include <cstdlib>
#include <unordered_map>
#include <armadillo>


namespace libwfa {

    class OmDescriptor {
    public:
        std::unordered_map<std::string, double> descriptor; //!< Store descriptor data using Python dict style

    private:
        const double om_tot; //!< Total omega
        const arma::mat &om_frag; //!< Omega matrix
        arma::mat om_norm; //!< normalized matrix
        size_t num_frag; //!< Number of fragments

    public:
        OmDescriptor(const double tot_om, const arma::mat &frag_om);

        virtual ~OmDescriptor() = default;

        double ret_desc(const std::string &desc);

        void ret_desc(const std::vector<std::string> &descs);

        void compute_desc(const std::string &desc);

        void compute_trans_met();

    }; //OmDescriptor

} //namespace

#endif //LIBWFA_OM_DESCRIPTOR_H
