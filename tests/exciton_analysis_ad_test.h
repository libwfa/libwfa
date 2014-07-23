#ifndef LIBWFA_EXCITON_ANALYSIS_AD_TEST_H
#define LIBWFA_EXCITON_ANALYSIS_AD_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/** \brief Tests for libwfa::exciton_analysis_ad class

    \ingroup libwfa_tests
 **/
class exciton_analysis_ad_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
private:
    void test_1();
};


}// end namespace libwfa


#endif  // LIBWFA_EXCITON_ANALYSIS_AD_TEST_H
