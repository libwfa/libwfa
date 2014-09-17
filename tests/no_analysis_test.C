#include <fstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/no_analysis.h>
#include "no_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void no_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
}

template<typename TestData>
void no_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "no_analysis_test::test_1()";

    try {

    std::ofstream of("no_analysis_test", std::ofstream::app);
    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = mat(nao, nmo);
    if (! data.aeqb()) c.beta() = mat(nao, nmo);
    read_ab_matrix(data, testname, "c", c);
        
    for (size_t istate = 1; istate <= data.nstates(); istate++) {

        ab_matrix dm(data.aeqb());
        dm.alpha() = mat(nao, nao);
        if (! data.aeqb()) dm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "dm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), dm);

        no_analysis na(s, c, dm);
        na.analyse(of);
        of << std::endl;
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
