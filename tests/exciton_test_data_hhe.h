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
    static arma::mat overlap() { return arma::mat(k_s,4,4); }
    static arma::mat x()  { return arma::mat(k_x, 4, 4); }
    static arma::mat y()  { return arma::mat(k_y, 4, 4); }
    static arma::mat z()  { return arma::mat(k_z, 4, 4); }
    static arma::mat xx() { return arma::mat(k_xx, 4, 4); }
    static arma::mat yy() { return arma::mat(k_yy, 4, 4); }
    static arma::mat zz() { return arma::mat(k_zz, 4, 4); }
    static arma::mat tdm() { return arma::mat(k_tdm, 4, 4); }
    static arma::vec rh() { return arma::vec(k_ref, 3); }
    static arma::vec re() { return arma::vec(k_ref + 3, 3); }
    static arma::vec rhre() { return arma::vec(k_ref + 6, 3); }
    static arma::vec rh2() { return arma::vec(k_ref + 9, 3); }
    static arma::vec re2() { return arma::vec(k_ref + 12, 3); }
};


}// end namespace libwfa

#endif  // LIBWFA_EXCITON_TEST_DATA_HHE_H
