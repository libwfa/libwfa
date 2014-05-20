#include "libwfa_suite.h"

namespace libwfa {


libwfa_suite::libwfa_suite() : libtest::test_suite("libwfa") {

    add_test("ab_matrix", m_utf_ab_matrix);
    add_test("ab_selector", m_utf_ab_selector);
    add_test("ab_vector", m_utf_ab_vector);
    add_test("ctnum_analysis", m_utf_ctnum_analysis);
    add_test("ctnum_export", m_utf_ctnum_export);
    add_test("export_densities_cube", m_utf_export_densities_cube);
    add_test("export_orbitals_cube", m_utf_export_orbitals_cube);
    add_test("export_orbitals_molden", m_utf_export_orbitals_molden);
    add_test("grid3d", m_utf_grid3d);
    add_test("pop_mulliken", m_utf_pop_mulliken);
    add_test("pop_print_default", m_utf_pop_print_default);
    add_test("selector", m_utf_selector);
    add_test("transformations_dm", m_utf_transformations_dm);
    add_test("version", m_utf_version);
}


} // namespace libwfa
