#include "libwfa_suite.h"

namespace libwfa {


libwfa_suite::libwfa_suite() : libtest::test_suite("libwfa") {

    add_test("ab_matrix", m_utf_ab_matrix);
    add_test("ab_vector", m_utf_ab_vector);
    add_test("pop_mulliken", m_utf_pop_mulliken);
    add_test("pop_print_default", m_utf_pop_print_default);
    add_test("transformations_dm", m_utf_transformations_dm);
    add_test("version", m_utf_version);
}


} // namespace libwfa
