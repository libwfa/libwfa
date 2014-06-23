#ifndef LIBWFA_EX_ANALYSE_AD_TEST_H
#define LIBWFA_EX_ANALYSE_AD_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/**	\brief Tests for libwfa::ex_analyse_ad class

	\ingroup libwfa_tests
 **/

class ex_analyse_ad_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);
private:
	/**\brief Tests all exciton tests.
	 **/
	void test_ex_total();
}; //end class


}// end namespace libwfa


#endif  // LIBWFA_EX_ANALYSE_AD_TEST_H
