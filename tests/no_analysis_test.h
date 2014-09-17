#ifndef LIBWFA_NO_ANALYSIS_TEST_H
#define LIBWFA_NO_ANALYSIS_TEST_H

#include "test_base.h"

namespace libwfa {

class test_base;


/** \brief Tests the libwfa::no_analysis class

    \ingroup libwfa_tests
 **/
class no_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);
    
private:
    template<typename TestData>
    void test_1() throw(libtest::test_exception);    
};


} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_TEST_H
