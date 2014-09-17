#ifndef LIBWFA_NDO_ANALYSIS_TEST_H
#define LIBWFA_NDO_ANALYSIS_TEST_H

#include "test_base.h"

namespace libwfa {


/** \brief Tests the libwfa::ndo_analysis class

    \ingroup libwfa_tests
 **/
class ndo_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);
    
private:
    void test_form_ad_1a() throw(libtest::test_exception);

    void test_form_ad_1b() throw(libtest::test_exception);

    template<typename TestData>
    void test_form_ad_2() throw(libtest::test_exception);

    template<typename TestData>
    void test_1() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_NDO_ANALYSIS_TEST_H
