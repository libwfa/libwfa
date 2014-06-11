#ifndef LIBWFA_ORBITAL_SELECTOR_TEST_H
#define LIBWFA_ORBITAL_SELECTOR_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::orbital_selector class

	\ingroup libmo_tests
 **/
class orbital_selector_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

private:
	void test_1();
    void test_2();
    void test_3();
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_SELECTOR_TEST_H
