#ifndef LIBWFA_GRID3D_TEST_H
#define LIBWFA_GRID3D_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::grid3d class

    \ingroup libmo_tests
 **/
class grid3d_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1a();
    void test_1b();
    void test_1c();
    void test_2a();
    void test_2b();
    void test_3();
    void test_4();
};


} // namespace libwfa

#endif // LIBWFA_GRID3D_TEST_H
