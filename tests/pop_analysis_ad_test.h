#ifndef LIBWFA_POP_ANALYSIS_AD_TEST_H
#define LIBWFA_POP_ANALYSIS_AD_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::pop_analysis_ad class

    \ingroup libwfa_tests
 **/
class pop_analysis_ad_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    template<typename TestData>
    void test_1() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_AD_TEST_H
