#include <iomanip>
#include <libwfa/analyses/pop_mulliken.h>
#include "pop_mulliken_test.h"
#include "test01_data.h"

namespace libwfa {

using namespace arma;

void pop_mulliken_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
}


void pop_mulliken_test::test_1() {

    static const char *testname = "pop_mulliken_test::test_1()";

    try {

    size_t na = 2, nb = 10, nb1 = 6;

    Col<size_t> b2c(nb, fill::zeros);
    for (size_t i = nb1; i < nb; i++) b2c(i) = 1;

    // Use the upper and lower triagonal of a random matrix to
    // form symmetric overlap and density matrices
    Mat<double> base = randu< Mat<double> >(nb, nb);
    Mat<double> ov = symmatu(base);
    Mat<double> dm = symmatl(base);

    Col<double> p, p_ref(na, fill::zeros);
    for (size_t i = 0; i < nb1; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++) {
            tmp += dm(i, j) * ov(i, j);
        }
        p_ref(0) -= tmp;
    }
    for (size_t i = nb1; i < nb; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++) {
            tmp += dm(i, j) * ov(i, j);
        }
        p_ref(1) -= tmp;
    }

    pop_mulliken(ov, b2c).perform(dm, p);

    if (p.size() != na) {
        fail_test(testname, __FILE__, __LINE__, "Length of population vector");
    }

    uvec x = find(abs(p - p_ref) > 1e-14, 1);
    if (x.size() != 0) {

        std::ostringstream oss;
        oss << "Population of atom " << x(0) << "(diff: " <<
                std::setprecision(6) << std::scientific <<
                p(x(0)) - p_ref(x(0)) << ")";
        fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


void pop_mulliken_test::test_2() {

    static const char *testname = "pop_mulliken_test::test_2()";

    try {

	size_t nao = test01_data::k_nao;
	size_t nmo = test01_data::k_nmo;
	test01_data data;
	Mat<double> s(nao, nao);
	ab_matrix dm(nao, nao, nao, nao);
	Col<size_t> b2p(nao);
	Col<double> pa_ref(test01_data::k_natoms), pa;

	data.read_matrix(testname, "s", s);
	data.read_ab_matrix(testname, "dm1", dm);
	data.read_matrix(testname, "b2at", b2p);
//	data.read_matrix(testname, "p1mulliken", pa_ref);

	pop_mulliken(s, b2p).perform(dm.alpha(), pa);

	if (pa.n_elem != pa_ref.n_elem) {
        fail_test(testname, __FILE__, __LINE__, "Length of population vector");
	}

//	uvec x = find(abs(p - p_ref) > 1e-14, 1);
//    if (x.size() != 0) {
//
//        std::ostringstream oss;
//        oss << "Population of atom " << x(0) << "(diff: " <<
//                std::setprecision(6) << std::scientific <<
//                p(x(0)) - p_ref(x(0)) << ")";
//        fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
//    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
