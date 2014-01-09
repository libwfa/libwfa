# 
# Module to find Armadillo C++ linear algebra library
# 
# Input variables:
# - ARMADILLO_PATH -- if this variable is defined, its value is used as the 
# 		first search path. 
#
# Output variables:
# - ARMADILLO_FOUND -- true or false, if the library has been found
# - ARMADILLO_INCLUDE_DIRS -- path(s) where to find the armadillo library
# - ARMADILLO_LIBRARIES -- the library to link against
# - ARMA_VERSION -- library version number 

set(ARMADILLO_FOUND FALSE)

#
#   1. Inspect ARMADILLO_PATH location
#
if (ARMADILLO_PATH) 
	find_library(ARMADILLO_LIBRARY NAMES armadillo 
		PATHS "${ARMADILLO_PATH}/lib" NO_DEFAULT_PATH)
	find_path(ARMADILLO_INCLUDE_DIR NAMES armadillo
		PATHS "${ARMADILLO_PATH}/include" NO_DEFAULT_PATH)  
endif(ARMADILLO_PATH)

#
#   2. Inspect standard locations
#
if (NOT ARMADILLO_INCLUDE_DIR) 
	find_library(ARMADILLO_LIBRARY NAMES armadillo)
	find_path(ARMADILLO_INCLUDE_DIR NAMES armadillo)
endif(NOT ARMADILLO_INCLUDE_DIR)

#
#   3. Extract version information from header file (if available)
#
set(ARMA_VERSION_HPP 
	"${ARMADILLO_INCLUDE_DIR}/armadillo_bits/arma_version.hpp")

if(EXISTS "${ARMA_VERSION_HPP}")
	file(STRINGS "${ARMA_VERSION_HPP}" ARMA_VERSION_CONTENTS 
		REGEX "^[ \t]*#define[ \t]+ARMA_VERSION_")
	string(REGEX REPLACE ".*ARMA_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" 
		ARMA_VERSION_MAJOR "${ARMA_VERSION_CONTENTS}")
	string(REGEX REPLACE ".*ARMA_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" 
		ARMA_VERSION_MINOR "${ARMA_VERSION_CONTENTS}")
	string(REGEX REPLACE ".*ARMA_VERSION_PATCH[ \t]+([0-9]+).*" "\\1" 
		ARMA_VERSION_PATCH "${ARMA_VERSION_CONTENTS}")
	set(ARMADILLO_VERSION 
		"${ARMA_VERSION_MAJOR}.${ARMA_VERSION_MINOR}.${ARMA_VERSION_PATCH}")
	set(ARMA_VERSION_HPP)
	set(ARMA_VERSION_MAJOR)
	set(ARMA_VERSION_MINOR)
	set(ARMA_VERSION_PATCH)
else() 
	set(ARMADILLO_VERSION "0.0.0")
endif()
	
message(STATUS "ARMADILLO_VERSION: ${ARMADILLO_VERSION}")

#
#   4. Set output variables
#
if(ARMADILLO_INCLUDE_DIR)
	get_filename_component(ARMADILLO_PATH "${ARMADILLO_INCLUDE_DIR}" PATH)
	set(ARMADILLO_INCLUDE_DIRS "${ARMADILLO_INCLUDE_DIR}")
	add_library(armadillo UNKNOWN IMPORTED)
	set_target_properties(armadillo PROPERTIES
        IMPORTED_LOCATION ${ARMADILLO_LIBRARY})
	set(ARMADILLO_LIBRARIES armadillo)	
	set(ARMADILLO_FOUND TRUE)
endif(ARMADILLO_INCLUDE_DIR)

#
#	5. Unset temporary variables
#
unset(ARMADILLO_LIBRARY)
unset(ARMADILLO_INCLUDE_DIR)
unset(ARMA_VERSION_HPP)
