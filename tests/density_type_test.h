#ifndef LIBWFA_DENSITY_TYPE_TEST_H
#define LIBWFA_DENSITY_TYPE_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::density_type class

	\ingroup libwfa_tests
 **/
class density_type_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_DENSITY_TYPE_TEST_H
