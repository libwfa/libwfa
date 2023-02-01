#ifndef LIBWFA_EXCITON_ANALYSIS_TEST_H
#define LIBWFA_EXCITON_ANALYSIS_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/**    \brief Tests for libwfa::exciton_analysis class

    \ingroup libwfa_tests
 **/
class exciton_analysis_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    /** \brief restricted case
     *
     */
    void test_1();

    /** \brief unrestricted case
     *
     */
    void test_2();
};


}// end namespace libwfa


#endif  // LIBWFA_EXCITON_ANALYSIS_TEST_H
