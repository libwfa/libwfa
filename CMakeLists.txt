cmake_minimum_required(VERSION 2.6)
project(LIBWFA)

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


if (${HAVE_LIBTEST})
	add_subdirectory(tests)
else (${HAVE_LIBTEST})
	Message(WARNING "libtest not found. Just compiling the bare library, but no tests.")
endif(${HAVE_LIBTEST})

