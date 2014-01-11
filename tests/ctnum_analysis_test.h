#ifndef LIBWFA_CTNUM_ANALYSIS_TEST_H
#define LIBWFA_CTNUM_ANALYSIS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::ctnum_analysis class

	\ingroup libwfa_tests
 **/
class ctnum_analysis_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

private:
	void test_1();

};


} // namespace libwfa

#endif // LIBWFA_CTNUM_ANALYSIS_TEST_H
