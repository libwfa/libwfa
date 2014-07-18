#ifndef LIBWFA_SANTO_ANALYSIS_TEST_H
#define LIBWFA_SANTO_ANALYSIS_TEST_H

#include "test_base.h"
#include <libwfa/analyses/santo_analysis.h>

namespace libwfa {

using namespace arma;
    
/** \brief Tests the libwfa::santo_analysis class

    \ingroup libwfa_tests
 **/
class santo_analysis_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    template<typename TestData>
    void test_1() throw(libtest::test_exception);

    void check(const ab_matrix &tdm, const ab_matrix &ui, const ab_matrix &vit,
               const ab_matrix &x, const mat &s, const char* testname) {
        
        check(tdm.alpha(), ui.alpha(), vit.alpha(), x.alpha(), s, testname);
        check(tdm.beta(),  ui.beta(),  vit.beta() , x.beta(), s, testname);
    }
    
    void check(const mat &tdm_x, const mat &ui_x, const mat &vit_x,
               const mat &x_x, const mat &s, const char* testname) {
        
        //Mat<double> tdm_chk = u_x.t() * tdm_x * v_x;
        
        if (accu(abs(ui_x.t() * ui_x - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "U not unitary.");
        if (accu(abs(vit_x * vit_x.t() - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "V not unitary.");
        //if (accu(abs(ev_chk_x % ev_chk_x - x_x) > 1e-12) != 0)
          //  fail_test(testname, __FILE__, __LINE__, "Bad transform.");
        
    }
};


} // namespace libwfa

#endif // LIBWFA_SANTO_ANALYSIS_TEST_H
