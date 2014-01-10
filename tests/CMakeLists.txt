set(SRC_HEADER
	ab_matrix_test.h
	ab_vector_test.h
	libwfa_suite.h
	pop_mulliken_test.h
	pop_print_default_test.h
    version_test.h
)

set(SRC
	ab_matrix_test.C
	ab_vector_test.C
    libwfa_tests.C
    libwfa_suite.C
	pop_mulliken_test.C
	pop_print_default_test.C
    version_test.C
)
add_executable(libwfa_tests ${SRC} ${SRC_HEADER})

include(${LIBTEST_DIR}/cmake/SetupLibTest.txt)

include_directories(${LIBWFA_DIR})

set_target_properties(libwfa_tests PROPERTIES
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

target_link_libraries(libwfa_tests wfa ${LIBWFA_DEP_LIBS} test)
