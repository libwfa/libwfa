#ifndef LIBWFA_VERSION_TEST_H
#define LIBWFA_VERSION_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/**	\brief Tests the libwfa::version class

	\ingroup libwfa_tests
 **/
class version_test : public libtest::unit_test {
public:
	virtual void perform() throw(libtest::test_exception);

};


} // namespace libwfa

#endif // LIBWFA_VERSION_TEST_H
