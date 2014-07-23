#ifndef LIBWFA_NDO_ANALYSIS_TEST_H
#define LIBWFA_NDO_ANALYSIS_TEST_H

#include "test_base.h"
#include <libwfa/analyses/ndo_analysis.h>

namespace libwfa {


/** \brief Tests the libwfa::ndo_analysis class

    \ingroup libwfa_tests
 **/
class ndo_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);
    
private:
    template<typename TestData>
    void test_1() throw(libtest::test_exception);    
};


} // namespace libwfa

#endif // LIBWFA_NDO_ANALYSIS_TEST_H
