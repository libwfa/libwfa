#include <iomanip>
#include <libwfa/analyses/pop_loewdin.h>
#include "pop_loewdin_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;

void pop_loewdin_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2<test01_data>();
    test_2<test02_data>();
}


void pop_loewdin_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "pop_loewdin_test::test_1()";

    try {

    size_t na = 2, nb = 10, nb1 = 6;

    uvec b2c(nb, fill::zeros);
    for (size_t i = nb1; i < nb; i++) b2c(i) = 1;

    // Use the upper and lower triagonal of a random matrix to
    // form symmetric overlap and density matrices
    mat base(nb, nb, fill::randu);
    mat dm = symmatl(base);
    mat ov, ov2;
    {
    mat u;
    vec e;
    eig_sym(e, u, symmatu(base));
    ov = u * diagmat(abs(e)) * u.t();
    ov2 = u * diagmat(sqrt(abs(e))) * u.t();
    }

    vec p, p_ref(na, fill::zeros);
    for (size_t i = 0; i < nb1; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++)
        for (size_t k = 0; k < nb; k++) {
            tmp += ov2(i, j) * dm(j, k) * ov2(i, k);
        }
        p_ref(0) -= tmp;
    }
    for (size_t i = nb1; i < nb; i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < nb; j++)
        for (size_t k = 0; k < nb; k++) {
            tmp += ov2(i, j) * dm(j, k) * ov2(i, k);
        }
        p_ref(1) -= tmp;
    }

    pop_loewdin(ov, b2c).perform(dm, p);

    if (p.size() != na) {
        fail_test(testname, __FILE__, __LINE__, "Length of population vector");
    }

    uvec x = find(abs(p - p_ref) > 1e-13, 1);
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
void pop_loewdin_test::test_2() throw(libtest::test_exception) {

    static const char *testname = "pop_loewdin_test::test_2()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;
    uvec b2p = data.bf2nuclei();

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    for (size_t i = 0; i <= data.nstates(); i++) {

        ab_matrix dm(data.aeqb());
        dm.alpha() = mat(nao, nao);
        if (! data.aeqb()) dm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "dm" << i;
        read_ab_matrix(data, testname, ssdm.str().c_str(), dm);

        //vec pa, pa_ref(TestData::k_natoms);
        vec pa(TestData::k_natoms);
        vec pa_ref = data.pop_loewdin(i, true);

        pop_loewdin(s, b2p).perform(dm.alpha(), pa);

        if (pa.n_elem != pa_ref.n_elem) {
            fail_test(testname, __FILE__, __LINE__,
                    "Length of population vector");
        }

        uvec x = find(abs(pa - pa_ref) > 1e-13, 1);
        if (x.size() != 0) {

            std::ostringstream oss;
            oss << "\n State " << i << ", population of atom " << x(0)
                    << "(diff: " << std::setprecision(6) << std::scientific
                    << pa(x(0)) - pa_ref(x(0)) << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
        }
        
        if (! data.aeqb()) {
            vec pb(TestData::k_natoms);
            vec pb_ref = data.pop_loewdin(i, false);
            pop_loewdin(s, b2p).perform(dm.beta(), pb);

            if (pb.n_elem != pb_ref.n_elem) {
                fail_test(testname, __FILE__, __LINE__,
                        "Length of population vector");
            }

            uvec xb = find(abs(pb - pb_ref) > 1e-13, 1);
            if (xb.size() != 0) {

                std::ostringstream oss;
                oss << "\n State " << i << ", population of atom " << xb(0)
                        << "(diff: " << std::setprecision(6) << std::scientific
                        << pb(xb(0)) - pb_ref(xb(0)) << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
            }            
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
