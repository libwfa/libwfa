#ifndef LIBWFA_CTNUM_EXPORT_TEST_H
#define LIBWFA_CTNUM_EXPORT_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::ctnum_export class

	\ingroup libwfa_tests
 **/
class ctnum_export_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

private:
	void test_1();
    void test_2();
    void test_exc();
};


} // namespace libwfa

#endif // LIBWFA_CTNUM_EXPORT_TEST_H
