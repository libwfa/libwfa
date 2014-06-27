#ifndef LIBWFA_TRANSFORMATIONS_DM_TEST_H
#define LIBWFA_TRANSFORMATIONS_DM_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::transformations_dm class

	\ingroup libmo_tests
 **/
class transformations_dm_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

private:
	void test_form_eh_1a();
    void test_form_eh_1b();
    void test_form_om_1a();
    void test_form_om_1b();
    void test_diagonalize_dm_1a();
    void test_diagonalize_dm_1b();
    void test_form_ad_1a();
    void test_form_ad_1b();
    void test_form_ad_2();
    void test_form_ad_3();
};


} // namespace libwfa

#endif // LIBWFA_TRANSFORMATIONS_DM_TEST_H
