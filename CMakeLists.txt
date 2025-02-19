



cmake_minimum_required(VERSION 3.14...3.24)

project(qfr
        LANGUAGES CXX
        DESCRIPTION "MQT QFR  - A library for Quantum Functionality Representation"
        )

include(thesis.cmake)

# enable organization of targets into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# configuration options
option(DEPLOY "Configure for deployment")
option(BINDINGS "Configure for building Python bindings")
option(COVERAGE "Configure for coverage report generation")
option(GENERATE_POSITION_INDEPENDENT_CODE "Generate position independent code")
option(BUILD_QFR_TESTS "Also build tests for QFR project")

if (DEFINED ENV{DEPLOY})
	set(DEPLOY $ENV{DEPLOY} CACHE BOOL "Use deployment configuration from environment" FORCE)
	message(STATUS "Setting deployment configuration to '${DEPLOY}' from environment")
endif ()

# set deployment specific options
if (DEPLOY)
	# set the macOS deployment target appropriately
	set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "" FORCE)
endif ()

# build type settings
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
	set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

macro(check_submodule_present MODULENAME)
	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${MODULENAME}/CMakeLists.txt")
		message(FATAL_ERROR "${MODULENAME} submodule not cloned properly. Please run `git submodule update --init --recursive` from the main project directory")
	endif()
endmacro()

check_submodule_present(dd_package)
check_submodule_present(json)
check_submodule_present(pybind11)
check_submodule_present(pybind11_json)
check_submodule_present(zx)

# add main library code
add_subdirectory(src)

# add test code
if(BUILD_QFR_TESTS)
	enable_testing()
	include(GoogleTest)
	add_subdirectory(test)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	add_subdirectory(apps)
endif()

# add Python binding code
if(BINDINGS)
	add_subdirectory(mqt/qfr)
endif()