#include "om_descriptor.h"

namespace libwfa {

    using namespace arma;

    OmDescriptor::OmDescriptor(const double &tot_om, const arma::mat &frag_om) : om_tot(tot_om), om_frag(frag_om) {
        om_norm = om_frag / om_tot;
        num_frag = om_frag.n_rows;
    }

    double OmDescriptor::ret_desc(std:: string desc) {
        if ( descriptor.find(desc) == descriptor.end() )
             compute_desc(desc);

        return descriptor[desc];

    }

    void OmDescriptor::compute_desc(std::string desc) {

        if (desc == "POSi") {

            colvec c1 = sum(om_norm, 1);
            colvec c2 = linspace<colvec>(1, num_frag, num_frag);
            descriptor[desc] = accu(c1 % c2);

        }

        else if (desc == "POSf") {

            rowvec c1 = sum(om_norm, 0);
            rowvec c2 = linspace<rowvec>(1, num_frag, num_frag);
            descriptor[desc] = accu(c1 % c2);

        }

        std::cout << desc << " ===== " <<  descriptor[desc] << std::endl;
    }

}; //namespace

