#ifndef LIBWFA_SA_NTO_ANALYSIS_TEST_H
#define LIBWFA_SA_NTO_ANALYSIS_TEST_H

#include <libwfa/analyses/sa_nto_analysis.h>
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

    void check(const ab_matrix &tdm, const ab_matrix &ui, const ab_matrix &vit,
            const ab_matrix &x, const mat &s, const char* testname);
    
    void check(const mat &tdm_x, const mat &ui_x, const mat &vit_x,
            const mat &x_x, const mat &s, const char* testname);
};


} // namespace libwfa

#endif // LIBWFA_SA_NTO_ANALYSIS_TEST_H
