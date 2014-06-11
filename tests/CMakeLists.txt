set(SRC_HEADER
    ab_matrix_test.h
    ab_orbital_selector_test.h
    ab_vector_test.h
    ctnum_analysis_test.h
    ctnum_export_test.h
    ctnumbers_test.h
    cube_writer_test.h
    export_cube_base_test.h
    grid3d_test.h
    libwfa_suite.h
    ndo_analysis_test.h
    no_analysis_test.h
    nto_analysis_test.h
    orbital_selector_test.h
    pop_analysis_ad_test.h
    pop_analysis_dm_test.h
    pop_mulliken_test.h
    pop_printer_default_test.h
    selector_test.h
    transformations_dm_test.h
    version_test.h
)

set(SRC
    ab_matrix_test.C
    ab_orbital_selector_test.C
    ab_vector_test.C
    ctnum_analysis_test.C
    ctnum_export_test.C
    ctnumbers_test.C
    cube_writer_test.C
    export_cube_base_test.C
    grid3d_test.C
    libwfa_tests.C
    libwfa_suite.C
    ndo_analysis_test.C
    no_analysis_test.C
    nto_analysis_test.C
    orbital_selector_test.C
    pop_analysis_ad_test.C
    pop_analysis_dm_test.C
    pop_mulliken_test.C
    pop_printer_default_test.C
    selector_test.C
    transformations_dm_test.C
    version_test.C
)

set(TARGETS libwfa_tests)
add_executable(libwfa_tests ${SRC} ${SRC_HEADER})

include(${LIBTEST_DIR}/cmake/SetupLibTest.txt)

#include_directories(${LIBWFA_DIR})

set_property(TARGET ${TARGETS} APPEND PROPERTY
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

foreach (TARGET ${TARGETS})
    target_link_libraries(${TARGET} wfa ${LIBWFA_DEP_LIBS} test)
endforeach(TARGET)