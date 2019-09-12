#ifndef LIBWFA_CTNUMBERS_TEST_H
#define LIBWFA_CTNUMBERS_TEST_H

#include "test_base.h"
#include <libwfa/analyses/ctnum_analysis.h>
#include <libwfa/analyses/ctnumbers.h>
#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::ctnumbers class

    \ingroup libwfa_tests
 **/
class ctnumbers_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    template<typename TestData>
    void test_1(std::string ctnum_type) throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_CTNUMBERS_TEST_H
