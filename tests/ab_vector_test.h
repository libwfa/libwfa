#ifndef LIBWFA_AB_VECTOR_TEST_H
#define LIBWFA_AB_VECTOR_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::ab_vector class

    \ingroup libwfa_tests
 **/
class ab_vector_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1a();
    void test_1b();
    void test_2();
    void test_3();
    void test_4a();
    void test_4b();
};


} // namespace libwfa

#endif // LIBWFA_AB_VECTOR_TEST_H
