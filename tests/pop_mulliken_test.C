#include <iomanip>
#include <libwfa/analyses/pop_mulliken.h>
#include "pop_mulliken_test.h"

namespace libwfa {

using namespace arma;

void pop_mulliken_test::perform() throw(libtest::test_exception) {

    test_1();
}


void pop_mulliken_test::test_1() {

    static const char *testname = "pop_mulliken_test::test_1()";

    try {

    size_t na = 2, nb = 10, nb1 = 6;

    std::vector<size_t> b2c(nb, 0);
    for (size_t i = nb1; i < nb; i++) b2c[i] = 1;
    std::vector<double> p0(na, 0.0);

    // Use the upper and lower triagonal of a random matrix to
    // form symmetric overlap and density matrices
    Mat<double> base = randu< Mat<double> >(nb, nb);
    Mat<double> ov = symmatu(base);
    Mat<double> dm = symmatl(base);

    std::vector<double> p, p_ref(na, 0.0);
    for (size_t i = 0; i < nb1; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++) {
            tmp += dm(i, j) * ov(i, j);
        }
        p_ref[0] -= tmp;
    }
    for (size_t i = nb1; i < nb; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++) {
            tmp += dm(i, j) * ov(i, j);
        }
        p_ref[1] -= tmp;
    }

    pop_mulliken(ov, b2c, p0).perform(dm, p);

    if (p.size() != na) {
        fail_test(testname, __FILE__, __LINE__, "Length of population vector");
    }
    for (size_t i = 0; i < na; i++) {
        if (fabs(p[i] - p_ref[i]) > 1e-14) {
            std::ostringstream oss;
            oss << "Population of atom " << i << "(diff: " <<
                    std::setprecision(6) << std::scientific <<
                    p[i] - p_ref[i] << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
