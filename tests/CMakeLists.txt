set(SRC_HEADER
    ab_matrix_test.h
    ab_vector_test.h
    compare_ref.h
    ctnum_analysis_test.h
    ctnum_export_test.h
    ctnumbers_test.h
    cube_writer_test.h
    density_type_test.h
    exciton_analysis_ad_test.h
    exciton_analysis_test.h
    exciton_test_data_hhe.h
    export_cube_base_test.h
    grid3d_test.h
    libwfa_exception_test.h
    libwfa_suite.h
    mom_builder_test.h
    ndo_analysis_test.h
    no_analysis_test.h
    nto_analysis_test.h
    sa_nto_analysis_test.h
    orbital_data_test.h
    orbital_selector_test.h
    test_orbital_printer.h
    orbital_type_test.h
    pop_analysis_ad_test.h
    pop_analysis_dm_test.h
    pop_analysis_tdm_test.h
    pop_data_test.h
    pop_loewdin_test.h
    pop_mulliken_test.h
    selector_test.h
    test_base.h
    test_data_base.h
    test00_data.h
    test01_data.h
    test02_data.h
    version_test.h
)

set(SRC
    ab_matrix_test.C
    ab_vector_test.C
    compare_ref.C
    ctnum_analysis_test.C
    ctnum_export_test.C
    ctnumbers_test.C
    cube_writer_test.C
    density_type_test.C
    exciton_analysis_ad_test.C
    exciton_analysis_test.C
    exciton_test_data_hhe.C
    export_cube_base_test.C
    grid3d_test.C
    libwfa_exception_test.C
    libwfa_tests.C
    libwfa_suite.C
    mom_builder_test.C
    ndo_analysis_test.C
    no_analysis_test.C
    nto_analysis_test.C
    sa_nto_analysis_test.C
    orbital_data_test.C
    orbital_selector_test.C
    test_orbital_printer.C
    orbital_type_test.C
    pop_analysis_ad_test.C
    pop_analysis_dm_test.C
    pop_analysis_tdm_test.C
    pop_data_test.C
    pop_loewdin_test.C
    pop_mulliken_test.C
    selector_test.C
    test_base.C
    test_data_base.C
    test00_data.C
    test01_data.C
    test02_data.C
    version_test.C
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
list (APPEND LIBWFA_DEP_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set(TARGETS libwfa_tests)
add_executable(libwfa_tests ${SRC} ${SRC_HEADER})

include(${LIBTEST_DIR}/cmake/SetupLibTest.txt)

#include_directories(${LIBWFA_DIR})

set_property(TARGET ${TARGETS} APPEND PROPERTY
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

foreach (TARGET ${TARGETS})
    target_link_libraries(${TARGET} wfa ${LIBWFA_DEP_LIBS} libtest)
endforeach(TARGET)
