#include <iomanip>
#include <libwfa/ctnum_analysis.h>
#include "ctnum_analysis_test.h"

namespace libwfa {

using namespace arma;

void ctnum_analysis_test::perform() throw(libtest::test_exception) {

    test_1();
}


void ctnum_analysis_test::test_1() {

    static const char *testname = "ctnum_analysis_test::test_1()";

    try {

    size_t na = 2, nb = 10, nb1 = 6;

    std::vector<size_t> b2c(nb, 0);
    for (size_t i = nb1; i < nb; i++) b2c[i] = 1;

    // Use the upper and lower triagonal of a random matrix to
    // form omega matrices
    Mat<double> base = randu< Mat<double> >(nb, nb);
    Mat<double> om_ao = symmatl(base);

    Mat<double> om_at, om_ref(na, na);
    om_ref.fill(0.0);
    for (size_t i = 0; i < nb1; i++) {
        for (size_t j = 0; j < nb1; j++) {
            om_ref(0, 0) += om_ao(i, j);
        }
    }
    for (size_t i = 0; i < nb1; i++) {
        for (size_t j = nb1; j < nb; j++) {
            om_ref(0, 1) += om_ao(i, j);
        }
    }
    for (size_t i = nb1; i < nb; i++) {
        for (size_t j = 0; j < nb1; j++) {
            om_ref(1, 0) += om_ao(i, j);
        }
    }
    for (size_t i = nb1; i < nb; i++) {
        for (size_t j = nb1; j < nb; j++) {
            om_ref(1, 1) += om_ao(i, j);
        }
    }

    ctnum_analysis(b2c).perform(om_ao, om_at);

    if (om_at.n_rows != na || om_at.n_cols != na) {
        fail_test(testname, __FILE__, __LINE__, "Size of CT number data");
    }
    for (size_t i = 0; i < na * na; i++) {
        if (fabs(om_at[i] - om_ref[i]) > 1e-14) {
            std::ostringstream oss;
            oss << "CT number of atom " << i / na << " and atom " << i % na <<
                    "(diff: " << std::setprecision(6) << std::scientific <<
                    om_at[i] - om_ref[i] << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
