#ifndef LIBWFA_CTNUM_ANALYSIS_TEST_H
#define LIBWFA_CTNUM_ANALYSIS_TEST_H

#include "test_base.h"

namespace libwfa {


/** \brief Tests the libwfa::ctnum_analysis class

    \ingroup libwfa_tests
 **/
class ctnum_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_form_om_1() throw(libtest::test_exception);

    void test_1() throw(libtest::test_exception);

};


} // namespace libwfa

#endif // LIBWFA_CTNUM_ANALYSIS_TEST_H
