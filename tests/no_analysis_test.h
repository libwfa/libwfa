#ifndef LIBWFA_NO_ANALYSIS_TEST_H
#define LIBWFA_NO_ANALYSIS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::no_analysis class

    \ingroup libwfa_tests
 **/
class no_analysis_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_TEST_H
