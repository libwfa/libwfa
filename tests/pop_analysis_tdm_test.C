#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/pop_analysis_tdm.h>
#include <libwfa/analyses/pop_mulliken.h>
#include "pop_analysis_tdm_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void pop_analysis_tdm_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
}


template<typename TestData>
void pop_analysis_tdm_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "pop_analysis_tdm_test::test_1()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;

    TestData data;
    uvec b2p = data.bf2nuclei();
    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    for (size_t i = 1; i <= data.nstates(); i++) {

        ab_matrix tdm(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        if (! data.aeqb()) tdm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "tdm" << i;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

        pop_mulliken pop(s, b2p);
        pop_data pd;

        pop_analysis_tdm(pop, tdm).perform(pd);

        if (pd.size() != 1) {
            fail_test(testname, __FILE__, __LINE__, "# data sets");
        }
        pop_data::iterator it = pd.begin();
        if (pd.name(it) != "Trans. (e)") {
            fail_test(testname, __FILE__, __LINE__, "Name of 1st data set");
        }
        double tot = accu(pd.data(it));
        if (fabs(tot) > 1e-12) {
            std::ostringstream oss;
            oss << "Non-zero total charge (" << tot << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
        }
        it++;
        if (it != pd.end()) {
            fail_test(testname, __FILE__, __LINE__, "it");
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
