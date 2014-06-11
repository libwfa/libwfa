#include "libwfa_suite.h"

namespace libwfa {


libwfa_suite::libwfa_suite() : libtest::test_suite("libwfa") {

    add_test("ab_matrix", m_utf_ab_matrix);
    add_test("ab_orbital_selector", m_utf_ab_orbital_selector);
    add_test("ab_vector", m_utf_ab_vector);
    add_test("ctnum_analysis", m_utf_ctnum_analysis);
    add_test("ctnum_export", m_utf_ctnum_export);
    add_test("ctnumbers", m_utf_ctnumbers);
    add_test("cube_writer", m_utf_cube_writer);
    add_test("export_cube_base", m_utf_export_cube_base);
    add_test("grid3d", m_utf_grid3d);
    add_test("ndo_analysis", m_utf_ndo_analysis);
    add_test("no_analysis", m_utf_no_analysis);
    add_test("nto_analysis", m_utf_nto_analysis);
    add_test("orbital_selector", m_utf_orbital_selector);
    add_test("pop_analysis_ad", m_utf_pop_analysis_ad);
    add_test("pop_analysis_dm", m_utf_pop_analysis_dm);
    add_test("pop_mulliken", m_utf_pop_mulliken);
    add_test("pop_printer_default", m_utf_pop_printer_default);
    add_test("selector", m_utf_selector);
    add_test("transformations_dm", m_utf_transformations_dm);
    add_test("version", m_utf_version);
}


} // namespace libwfa
