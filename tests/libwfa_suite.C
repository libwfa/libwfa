#include "libwfa_suite.h"

namespace libwfa {


libwfa_suite::libwfa_suite() : libtest::test_suite("libwfa") {

    add_test("pop_mulliken", m_utf_pop_mulliken);
    add_test("pop_print_default", m_utf_pop_print_default);
//    add_test("version", m_utf_version);
}


} // namespace libwfa
