cmake_minimum_required(VERSION 2.6)
project(LIBWFA)

get_filename_component(LIBTEST_DIR ${LIBWFA_SOURCE_DIR}/../libtest ABSOLUTE)
set(LIBWFA_DIR ${LIBWFA_SOURCE_DIR})

include(${LIBWFA_DIR}/cmake/DepsLibWFA.txt)

include_directories(${LIBWFA_DIR})

add_subdirectory(libwfa)
add_subdirectory(tests)

