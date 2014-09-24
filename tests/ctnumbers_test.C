#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include "ctnumbers_test.h"
#include "test00_data.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void ctnumbers_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
}


template<typename TestData>
void ctnumbers_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ctnumbers_test::test_1()";

    try {
        
    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    size_t na  = TestData::k_natoms;
    TestData data;
    uvec b2p = data.bf2nuclei();

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ctnum_analysis cta(s, b2p);

    for (size_t istate = 1; istate <= data.nstates(); istate++) {

        ab_matrix tdm(data.aeqb()), om_ref(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        om_ref.alpha() = mat(na, na);
        if (! data.aeqb()) {
            tdm.beta() = mat(nao, nao);
            om_ref.beta() = mat(na, na);
        }

        std::ostringstream ssdm; ssdm << "tdm" << istate;
        std::ostringstream ssom; ssom << "om" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);
        read_ab_matrix(data, testname, ssom.str().c_str(), om_ref);

        ctnumbers ctnum(cta, tdm);
        const ab_matrix &om = ctnum.omega();
        if (om.is_alpha_eq_beta() != data.aeqb())
            fail_test(testname, __FILE__, __LINE__, "Alpha == beta.");
        if (om.alpha().n_rows != na || om.alpha().n_cols != na)
            fail_test(testname, __FILE__, __LINE__, "Size of CT number data");
        if (accu(abs(om_ref.alpha() - om.alpha())) > 1e-12)
            fail_test(testname, __FILE__, __LINE__, "CT number (alpha)");
        if (data.aeqb()) {
            if (om.beta().n_rows != na || om.beta().n_cols != na)
                fail_test(testname, __FILE__, __LINE__, "Size of CT number data");
            if (accu(abs(om_ref.beta() - om.beta())) > 1e-12)
                fail_test(testname, __FILE__, __LINE__, "CT number (alpha)");
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
    
}

} // namespace libwfa
