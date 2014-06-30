#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/pop_analysis_ad.h>
#include <libwfa/analyses/pop_mulliken.h>
#include "pop_analysis_ad_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void pop_analysis_ad_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
}


template<typename TestData>
void pop_analysis_ad_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "pop_analysis_ad_test::test_1()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;

    TestData data;
    Col<size_t> b2p = data.bf2nuclei();
    Mat<double> s(nao, nao);
    data.read_matrix(testname, "s", s);

    for (size_t i = 0; i < data.nstates(); i++) {

        ab_matrix at(data.aeqb()), de(data.aeqb());
        at.alpha() = Mat<double>(nao, nao);
        de.alpha() = Mat<double>(nao, nao);
        if (! data.aeqb()) {
            at.beta() = Mat<double>(nao, nao);
            de.beta() = Mat<double>(nao, nao);
        }

        std::ostringstream ssat; ssat << "atdm" << i + 1;
        std::ostringstream ssde; ssde << "dedm" << i + 1;
        data.read_ab_matrix(testname, ssat.str().c_str(), at);
        data.read_ab_matrix(testname, ssde.str().c_str(), de);

        pop_mulliken pop(s, b2p);
        pop_data pd;
        pop_analysis_ad(pop, at, de).perform(pd);

        if (data.aeqb()) {
            if (pd.size() != 3) {
                fail_test(testname, __FILE__, __LINE__, "# data sets");
            }
            pop_data::iterator it = pd.begin();
            if (pd.name(it) != "h+") {
                fail_test(testname, __FILE__, __LINE__, "Name of 1st data set");
            }
            it++;
            if (pd.name(it) != "e-") {
                fail_test(testname, __FILE__, __LINE__, "Name of 2nd data set");
            }
            it++;
            if (pd.name(it) != "Del q") {
                fail_test(testname, __FILE__, __LINE__, "Name of 1st data set");
            }
            double tot = accu(pd.data(it));
            if (fabs(tot) > 1e-5) {
                std::ostringstream oss;
                oss << "Non-zero total (" << tot << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
            }
            it++;
            if (it != pd.end()) {
                fail_test(testname, __FILE__, __LINE__, "it");
            }
        }
        else {
            if (pd.size() != 4) {
                fail_test(testname, __FILE__, __LINE__, "# data sets");
            }
            pop_data::iterator it = pd.begin();
            if (pd.name(it) != "h+ (alpha)") {
                fail_test(testname, __FILE__, __LINE__, "Name of 1st data set");
            }
            it++;
            if (pd.name(it) != "h+ (beta)") {
                fail_test(testname, __FILE__, __LINE__, "Name of 2nd data set");
            }
            it++;
            if (pd.name(it) != "e- (alpha)") {
                fail_test(testname, __FILE__, __LINE__, "Name of 3rd data set");
            }
            it++;
            if (pd.name(it) != "e- (beta)") {
                fail_test(testname, __FILE__, __LINE__, "Name of 4th data set");
            }
            it++;
            if (it != pd.end()) {
                fail_test(testname, __FILE__, __LINE__, "it");
            }
        }
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
