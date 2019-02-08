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
        arma::mat om_norm; //!< norm matrix
        size_t num_frag; //!< Number of fragments

        public:
        OmDescriptor(const double &tot_om, const arma::mat &frag_om);

        virtual ~OmDescriptor() = default;

        double ret_desc(std:: string desc);

        void compute_desc(std:: string desc);

    }; //OmDescriptor

} //namespace

#endif //LIBWFA_OM_DESCRIPTOR_H
