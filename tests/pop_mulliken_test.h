#ifndef LIBWFA_POP_MULLIKEN_TEST_H
#define LIBWFA_POP_MULLIKEN_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::pop_mulliken class

    \ingroup libwfa_tests
 **/
class pop_mulliken_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1() throw(libtest::test_exception);

    template<typename TestData>
    void test_2() throw(libtest::test_exception);

};


} // namespace libwfa

#endif // LIBWFA_POP_MULLIKEN_TEST_H
