#ifndef LIBWFA_POP_ANALYSIS_DM_TEST_H
#define LIBWFA_POP_ANALYSIS_DM_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::pop_analysis_dm class

	\ingroup libwfa_tests
 **/
class pop_analysis_dm_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_DM_TEST_H
