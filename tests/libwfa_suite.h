#ifndef LIBWFA_LIBWFA_SUITE_H
#define LIBWFA_LIBWFA_SUITE_H

#include <libtest/test_suite.h>
#include "ab_matrix_test.h"
#include "ab_selector_test.h"
#include "ab_vector_test.h"
#include "ctnum_analysis_test.h"
#include "ctnum_export_test.h"
#include "export_densities_cube_test.h"
#include "export_orbitals_cube_test.h"
#include "export_orbitals_molden_test.h"
#include "grid3d_test.h"
#include "pop_mulliken_test.h"
#include "pop_printer_default_test.h"
#include "selector_test.h"
#include "transformations_dm_test.h"
#include "version_test.h"

using libtest::unit_test_factory;

/** \defgroup libwfa_tests Tests
    \ingroup libwfa
 **/

namespace libwfa {


/** \brief Test suite for the wfa library (libwfa)

    This suite runs the following tests:
    - \c ab_matrix_test
    - \c ab_selector_test
    - \c ab_vector_test
    - \c ctnum_analysis_test
    - \c ctnum_export_test
    - \c export_densities_cube_test
    - \c export_orbitals_cube_test
    - \c export_orbitals_molden_test
    - \c grid3d_test
    - \c pop_mulliken_test
    - \c pop_printer_default_test
    - \c selector_test
    - \c transformations_dm_test
    - \c version_test

    \ingroup libwfa_tests
 **/
class libwfa_suite: public libtest::test_suite {
private:
    unit_test_factory<ab_matrix_test> m_utf_ab_matrix;
    unit_test_factory<ab_selector_test> m_utf_ab_selector;
    unit_test_factory<ab_vector_test> m_utf_ab_vector;
    unit_test_factory<ctnum_analysis_test> m_utf_ctnum_analysis;
    unit_test_factory<ctnum_export_test> m_utf_ctnum_export;
    unit_test_factory<export_densities_cube_test> m_utf_export_densities_cube;
    unit_test_factory<export_orbitals_cube_test> m_utf_export_orbitals_cube;
    unit_test_factory<export_orbitals_molden_test> m_utf_export_orbitals_molden;
    unit_test_factory<pop_mulliken_test> m_utf_grid3d;
    unit_test_factory<pop_mulliken_test> m_utf_pop_mulliken;
    unit_test_factory<pop_printer_default_test> m_utf_pop_printer_default;
    unit_test_factory<selector_test> m_utf_selector;
    unit_test_factory<transformations_dm_test> m_utf_transformations_dm;
    unit_test_factory<version_test> m_utf_version;

public:
    //! Creates the suite
    libwfa_suite();

};


} // namespace libwfa

#endif // LIBWFA_LIBWFA_SUITE_H
