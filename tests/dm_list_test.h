#ifndef LIBWFA_DM_LIST_TEST_H
#define LIBWFA_DM_LIST_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::dm_list class

	\ingroup libmo_tests
 **/
class dm_list_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
    void test_2();
    void test_3();
    void test_4();
};


} // namespace libwfa

#endif // LIBWFA_DM_LIST_TEST_H
