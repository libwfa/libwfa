#ifndef LIBWFA_EXPORT_ORBITALS_MOLDEN_TEST_H
#define LIBWFA_EXPORT_ORBITALS_MOLDEN_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::export_orbitals_molden class

	\ingroup libmo_tests
 **/
class export_orbitals_molden_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
    void test_2();
    void test_3();
};


} // namespace libwfa

#endif // LIBWFA_EXPORT_ORBITALS_MOLDEN_TEST_H
