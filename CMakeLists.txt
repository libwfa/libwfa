cmake_minimum_required(VERSION 2.6)
project(LIBWFA)

option(MOLCAS_LIB  "Compile library for Molcas" OFF)
option(MOLCAS_EXE  "Compile standalone Molcas executable" OFF)
option(ARMA_HEADER "Directly use armadillo headers" OFF)

if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

if(MOLCAS_LIB)
    set(ARMA_HEADER ON)
endif(MOLCAS_LIB)

get_filename_component(LIBTEST_DIR ${LIBWFA_SOURCE_DIR}/../libtest ABSOLUTE)
if (IS_DIRECTORY ${LIBTEST_DIR})
	set(HAVE_LIBTEST TRUE)
else (IS_DIRECTORY ${LIBTEST_DIR})
	set(HAVE_LIBTEST FALSE)
endif (IS_DIRECTORY ${LIBTEST_DIR})

set(LIBWFA_DIR ${LIBWFA_SOURCE_DIR})

include(${LIBWFA_DIR}/cmake/DepsLibWFA.txt)

include_directories(${LIBWFA_DIR})

add_subdirectory(libwfa)


if (HAVE_LIBTEST)
	add_subdirectory(tests)
else (HAVE_LIBTEST)
	Message(WARNING "libtest not found. Just compiling the bare library, but no tests.")
endif(HAVE_LIBTEST)
