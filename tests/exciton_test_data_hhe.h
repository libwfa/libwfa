#ifndef LIBWFA_EXCITON_TEST_DATA_HHE_H
#define LIBWFA_EXCITON_TEST_DATA_HHE_H

#include <armadillo>

namespace libwfa {


/** \brief Test data for H-He test case
 **/
class exciton_test_data_hhe {
private:
    static const double k_s[16];
    static const double k_x[16];
    static const double k_y[16];
    static const double k_z[16];
    static const double k_xx[16];
    static const double k_yy[16];
    static const double k_zz[16];
    static const double k_tdm[16];
    static const double k_ref[15];

public:
    static arma::Mat<double> overlap() { return arma::Mat<double>(k_s,4,4); }
    static arma::Mat<double> x()  { return arma::Mat<double>(k_x, 4, 4); }
    static arma::Mat<double> y()  { return arma::Mat<double>(k_y, 4, 4); }
    static arma::Mat<double> z()  { return arma::Mat<double>(k_z, 4, 4); }
    static arma::Mat<double> xx() { return arma::Mat<double>(k_xx, 4, 4); }
    static arma::Mat<double> yy() { return arma::Mat<double>(k_yy, 4, 4); }
    static arma::Mat<double> zz() { return arma::Mat<double>(k_zz, 4, 4); }
    static arma::Mat<double> tdm() { return arma::Mat<double>(k_tdm, 4, 4); }
    static arma::Col<double> rh() { return arma::Col<double>(k_ref, 3); }
    static arma::Col<double> re() { return arma::Col<double>(k_ref + 3, 3); }
    static arma::Col<double> rhre() { return arma::Col<double>(k_ref + 6, 3); }
    static arma::Col<double> rh2() { return arma::Col<double>(k_ref + 9, 3); }
    static arma::Col<double> re2() { return arma::Col<double>(k_ref + 12, 3); }
};


}// end namespace libwfa

#endif  // LIBWFA_EXCITON_TEST_DATA_HHE_H
