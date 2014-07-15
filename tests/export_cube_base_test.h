#ifndef LIBWFA_EXPORT_CUBE_BASE_TEST_H
#define LIBWFA_EXPORT_CUBE_BASE_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::grid3d class

    \ingroup libwfa_tests
 **/
class export_cube_base_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_CUBE_BASE_TEST_H
