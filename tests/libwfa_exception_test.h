#ifndef LIBWFA_LIBWFA_EXCEPTION_TEST_H
#define LIBWFA_LIBWFA_EXCEPTION_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::libwfa_exception class

    \ingroup libwfa_tests
 **/
class libwfa_exception_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_LIBWFA_EXCEPTION_TEST_H
