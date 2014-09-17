#ifndef LIBWFA_MOM_BUILDER_TEST_H
#define LIBWFA_MOM_BUILDER_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::mom_builder class

    \ingroup libwfa_tests
 **/
class mom_builder_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
};


} // namespace libwfa

#endif // LIBWFA_MOM_BUILDER_TEST_H
