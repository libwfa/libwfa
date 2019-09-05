#ifndef LIBWFA_SOC
#define LIBWFA_SOC

#include <armadillo>

struct h_so {
    arma::Mat<std::complex<double> > lm1_sp1;
    arma::Mat<std::complex<double> > l0_s0;
    arma::Mat<std::complex<double> > lp1_sm1;
};

struct soc {
    std::complex<double> lm1_sp1;
    std::complex<double> l0_s0;
    std::complex<double> lp1_sm1;
};

#endif // LIBWFA_SOC
