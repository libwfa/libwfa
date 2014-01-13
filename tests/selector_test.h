#ifndef LIBWFA_SELECTOR_TEST_H
#define LIBWFA_SELECTOR_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::selector class

	\ingroup libmo_tests
 **/
class selector_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

private:
	void test_1();
    void test_2();
    void test_3();
    void test_4a();
    void test_4b();
    void test_5();
    void test_6();
    void test_7();
};


} // namespace libwfa

#endif // LIBWFA_SELECTOR_TEST_H
