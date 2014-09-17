#ifndef LIBWFA_SA_NTO_ANALYSIS_TEST_H
#define LIBWFA_SA_NTO_ANALYSIS_TEST_H

#include "test_base.h"

namespace libwfa {

using namespace arma;
    
/** \brief Tests the libwfa::sa_nto_analysis class

    \ingroup libwfa_tests
 **/
class sa_nto_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    template<typename TestData>
    void test_1() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_SA_NTO_ANALYSIS_TEST_H
