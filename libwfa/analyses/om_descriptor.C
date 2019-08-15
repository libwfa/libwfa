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


#include <exception>
#include "om_descriptor.h"

namespace libwfa {

    using namespace arma;

    OmDescriptor::OmDescriptor(const std::vector<double> &tot_om, const arma::mat &frag_om) : om_tot(tot_om), om_frag(frag_om) {
        om_norm = om_frag / om_tot[0];
        num_frag = om_frag.n_rows;
    }

    double OmDescriptor::ret_desc(const std::string &desc) {
        if ( descriptor.find(desc) == descriptor.end() ) {
            compute_desc(desc);
        }

        return descriptor[desc];

    }

    void OmDescriptor::ret_desc(const std::vector<std::string> &descs) {

        for (auto const& desc: descs) {
            if ( descriptor.find(desc) == descriptor.end() ) {
                compute_desc(desc);
            }
        }

    }

    void OmDescriptor::compute_desc(const std::string &desc) {

        if (desc == "Om") {

            descriptor[desc] = om_tot[0] + om_tot[1];

        }

        else if (desc == "POSi") {

            colvec c1 = sum(om_norm, 1);
            colvec c2 = linspace<colvec>(1, num_frag, num_frag);
            descriptor[desc] = accu(c1 % c2);

        }

        else if (desc == "POSf") {

            rowvec c1 = sum(om_norm, 0);
            rowvec c2 = linspace<rowvec>(1, num_frag, num_frag);
            descriptor[desc] = accu(c1 % c2);

        }

        else if (desc == "POS") {

            descriptor[desc] = 0.5 * (ret_desc("POSi") + ret_desc("POSf"));

        }

        else if (desc == "CT") {

            descriptor[desc] = 0.0;
            for ( size_t i = 0; i < num_frag; i++ ) {
                for ( size_t j = i + 1; j < num_frag; j++ ) {
                    descriptor[desc] += om_norm(i, j) + om_norm(j, i);
                }
            }
        }

        else if (desc == "CT2") {

            descriptor[desc] = 0.0;
            for ( size_t i = 0; i < num_frag; i++ ) {
                for ( size_t j = i + 2; j < num_frag; j++ ) {
                    descriptor[desc] += om_norm(i, j) + om_norm(j, i);
                }
            }
        }

        else if (desc == "CTnt") {

            descriptor[desc] = ret_desc("POSf") - ret_desc("POSi");

        }

        else if (desc == "PRi") {

            descriptor[desc] = 1. / accu( square( sum(om_norm, 1) ) );

        }

        else if (desc == "PRf" || desc == "EEDL") {

            descriptor[desc] = 1. / accu( square( sum(om_norm, 0) ) );

        }

        else if (desc == "PR") {

            descriptor[desc] = 0.5 * (ret_desc("PRi") + ret_desc("PRf"));

        }

        else if (desc == "PRh") {

            descriptor[desc] = 2. / (1.0 / ret_desc("PRi") + 1.0 / ret_desc("PRf"));

        }

        else if (desc == "DEL") {

            // method 1:: concise but more expensive
            descriptor[desc] = 1. / accu( square( sum( (om_norm + om_norm.t()) / 2., 1) ) ); // note om_norm.t() makes a copy of data

            // method 2:: looks ugly but efficient. Use method 2 if speed is a concern
            /*
            descriptor[desc] = 0.;
            for ( size_t i = 0; i < num_frag; i++ ) {
                double tmp = 0.0;

                for ( size_t j = 0; j < num_frag; j++ ) {
                    tmp += (om_norm(i, j) + om_norm(j, i)) / 2.;
                }

                descriptor[desc] += tmp * tmp;
            }
            descriptor[desc] = 1./ descriptor[desc]; */

        }

        else if (desc == "COH") {

            descriptor[desc] = 1. / accu( sum(square(om_norm), 1) ) / ret_desc("PR");

        }

        else if (desc == "COHh") {

            descriptor[desc] = 1. / accu( sum( square(om_norm), 1 ) ) /  ret_desc("PRh");

        }

        else if (desc == "MC" || desc == "LC" || desc == "MLCT" || desc == "LMCT" ||
                desc == "LLCT" || desc == "SIEL") {

            compute_trans_met();

        }

        else {

            throw("Unknown descriptor!");
        }

        //std::cout << "Descriptor " << desc << " = " <<  descriptor[desc] << std::endl;
    }

    void OmDescriptor::compute_trans_met() {

        descriptor["MC"] = om_norm(0, 0);

        descriptor["LC"]   = 0.;
        descriptor["MLCT"] = 0.;
        descriptor["LMCT"] = 0.;
        descriptor["LLCT"] = 0.;

        for ( size_t i = 1; i < num_frag; i++ ) {

            descriptor["LC"]   += om_norm(i, i);
            descriptor["MLCT"] += om_norm(0, i);
            descriptor["LMCT"] += om_norm(i, 0);

            for ( size_t j = i + 1; j < num_frag; j++ ) {

                descriptor["LLCT"] += om_norm(i, j) + om_norm(j, i);

            }

        }

        rowvec epop = sum(om_frag, 0);
        descriptor["SIEL"] = -1 * epop(1) + 0.5 * sum( epop.tail(epop.n_elem - 2) );

    }

}; //namespace
