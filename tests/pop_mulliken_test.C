#include <iomanip>
#include <libwfa/analyses/pop_mulliken.h>
#include "pop_mulliken_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;

void pop_mulliken_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2<test01_data>();
    test_2<test02_data>();
}


void pop_mulliken_test::test_1() throw(libtest::test_exception) {

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


template<typename TestData>
void pop_mulliken_test::test_2() throw(libtest::test_exception) {

    static const char *testname = "pop_mulliken_test::test_2()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;
    Col<size_t> b2p = data.bf2nuclei();

    Mat<double> s(nao, nao);
    read_matrix(data, testname, "s", s);

    for (size_t i = 0; i <= data.nstates(); i++) {

        ab_matrix dm(data.aeqb());
        dm.alpha() = Mat<double>(nao, nao);
        if (! data.aeqb()) dm.beta() = Mat<double>(nao, nao);

        std::ostringstream ssdm; ssdm << "dm" << i;
        read_ab_matrix(data, testname, ssdm.str().c_str(), dm);

        //Col<double> pa, pa_ref(TestData::k_natoms);
        Col<double> pa(TestData::k_natoms);
        Col<double> pa_ref = data.popref(i, true);

        pop_mulliken(s, b2p).perform(dm.alpha(), pa);

        if (pa.n_elem != pa_ref.n_elem) {
            fail_test(testname, __FILE__, __LINE__,
                    "Length of population vector");
        }

        uvec x = find(abs(pa - pa_ref) > 1e-14, 1);
        if (x.size() != 0) {

            std::ostringstream oss;
            oss << "\n State " << i << " , population of atom " << x(0) << "(diff: " <<
                    std::setprecision(6) << std::scientific <<
                    pa(x(0)) - pa_ref(x(0)) << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
        }
        
        if (! data.aeqb()) {
            Col<double> pb(TestData::k_natoms);
            Col<double> pb_ref = data.popref(i, false);
            pop_mulliken(s, b2p).perform(dm.beta(), pb);

            if (pb.n_elem != pb_ref.n_elem) {
                fail_test(testname, __FILE__, __LINE__,
                        "Length of population vector");
            }

            uvec xb = find(abs(pb - pb_ref) > 1e-14, 1);
            if (xb.size() != 0) {

                std::ostringstream oss;
                oss << "\n State " << i << " , population of atom " << xb(0) << "(diff: " <<
                        std::setprecision(6) << std::scientific <<
                        pb(xb(0)) - pb_ref(xb(0)) << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
            }            
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
