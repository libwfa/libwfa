#ifndef LIBWFA_CTNUMBERS_TEST_H
#define LIBWFA_CTNUMBERS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::ctnumbers class

    \ingroup libwfa_tests
 **/
class ctnumbers_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_TEST_H
