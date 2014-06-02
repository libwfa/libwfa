#ifndef EX_ANA_TEST_H_
#define EX_ANA_TEST_H_

#include <libtest/unit_test.h>

namespace libwfa {

/**	\brief Tests for libwfa::ex_analyse class

	\ingroup libmo_tests
 **/

class ex_ana_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);
private:
	/**\brief Tests all exciton tests and prints the results.
	 **/
	void test_ex_total();
}; //end class


}// end namespace libwfa


#endif  //EX_ANA_TEST_H_
