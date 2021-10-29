#ifndef LIBWFA_NON_ORTH_TEST_H
#define LIBWFA_NON_ORTH_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/**    \brief Tests for 1TDM from non-orthogonal orbitals

    \ingroup libwfa_tests
 **/
class non_orth_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();

    //void test_2();
};


}// end namespace libwfa


#endif  // LIBWFA_NON_ORTH_TEST_H
