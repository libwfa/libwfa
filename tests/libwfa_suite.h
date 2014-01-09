#ifndef LIBWFA_LIBWFA_SUITE_H
#define LIBWFA_LIBWFA_SUITE_H

#include <libtest/test_suite.h>
#include "pop_mulliken_test.h"
#include "pop_print_default_test.h"
//#include "version_test.h"

using libtest::unit_test_factory;

/** \defgroup libwfa_tests Tests
    \ingroup libwfa
 **/

namespace libwfa {


/** \brief Test suite for the wfaity library (libwfa)

    This suite runs the following tests:
    - \c pop_mulliken_test
    - \c pop_print_default_test

    \ingroup libwfa_tests
 **/
class libwfa_suite: public libtest::test_suite {
private:
    unit_test_factory<pop_mulliken_test> m_utf_pop_mulliken;
    unit_test_factory<pop_print_default_test> m_utf_pop_print_default;
//    unit_test_factory<version_test> m_utf_version;

public:
    //! Creates the suite
    libwfa_suite();

};


} // namespace libwfa

#endif // LIBWFA_LIBWFA_SUITE_H
