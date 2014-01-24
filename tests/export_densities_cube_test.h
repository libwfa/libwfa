#ifndef LIBWFA_EXPORT_DENSITIES_CUBE_TEST_H
#define LIBWFA_EXPORT_DENSITIES_CUBE_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::export_densities_cube class

	\ingroup libmo_tests
 **/
class export_densities_cube_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
    void test_2();
    void test_3();
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_DENSITIES_CUBE_TEST_H
