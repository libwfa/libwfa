#ifndef LIBWFA_AB_MATRIX_TEST_H
#define LIBWFA_AB_MATRIX_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::ab_matrix class

    \ingroup libwfa_tests
 **/
class ab_matrix_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1a();
    void test_1b();
    void test_2a();
    void test_2b();
    void test_3a();
    void test_3b();
    void test_4a();
    void test_4b();
};


} // namespace libwfa

#endif // LIBWFA_AB_MATRIX_TEST_H
