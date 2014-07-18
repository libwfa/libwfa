#ifndef LIBWFA_SANTO_ANALYSIS_TEST_H
#define LIBWFA_SANTO_ANALYSIS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::santo_analysis class

    \ingroup libwfa_tests
 **/
class santo_analysis_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_SANTO_ANALYSIS_TEST_H
