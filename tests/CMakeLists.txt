set(SRC
    libwfa_tests.C
    libwfa_suite.C
#    version_test.C
)
add_executable(libwfa_tests ${SRC})

include(${LIBTEST_DIR}/cmake/SetupLibTest.txt)

include_directories(${LIBWFA_DIR})

set_target_properties(libwfa_tests PROPERTIES
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

target_link_libraries(libwfa_tests wfa ${LIBWFA_DEP_LIBS} test)
