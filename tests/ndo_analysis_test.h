#ifndef LIBWFA_NDO_ANALYSIS_TEST_H
#define LIBWFA_NDO_ANALYSIS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::ndo_analysis class

    \ingroup libwfa_tests
 **/
class ndo_analysis_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_NDO_ANALYSIS_TEST_H
