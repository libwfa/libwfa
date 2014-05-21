set(SRC_HEADER
    ab_matrix_test.h
    ab_vector_test.h
    ab_selector_test.h
    ctnum_analysis_test.h
    ctnum_export_test.h
    libwfa_suite.h
    pop_mulliken_test.h
    pop_printer_default_test.h
    selector_test.h
    transformations_dm_test.h
    version_test.h
)

set(SRC
    ab_matrix_test.C
    ab_vector_test.C
    ab_selector_test.C
    ctnum_analysis_test.C
    ctnum_export_test.C
    libwfa_tests.C
    libwfa_suite.C
    pop_mulliken_test.C
    pop_printer_default_test.C
    selector_test.C
    transformations_dm_test.C
    version_test.C
)

add_executable(libwfa_tests ${SRC} ${SRC_HEADER})

include(${LIBTEST_DIR}/cmake/SetupLibTest.txt)

include_directories(${LIBWFA_DIR})

set_target_properties(libwfa_tests PROPERTIES
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

target_link_libraries(libwfa_tests wfa ${LIBWFA_DEP_LIBS} test)
