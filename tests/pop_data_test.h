#ifndef LIBWFA_POP_DATA_TEST_H
#define LIBWFA_POP_DATA_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::pop_data class

    \ingroup libwfa_tests
 **/
class pop_data_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
    void test_2();
    void test_exc();
};


} // namespace libwfa

#endif // LIBWFA_POP_DATA_TEST_H
