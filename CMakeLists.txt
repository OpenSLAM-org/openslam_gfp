# This is the main CMake configuration file, you should always do "cmake ." from here and not from subdirs

# Needed by Cmake versione 2.6 or later (Run "cmake --help-policy CMP0003" for more information.)"
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 OLD)
ENDIF(COMMAND CMAKE_POLICY)

SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};./build_tools")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
 
# Some project related variables
PROJECT(GFP CXX C)

SET(CMAKE_VERBOSE_MAKEFILE 0)
SET(VERBOSE 1)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${GFP_SOURCE_DIR}/src")
SET(EXECUTABLE_OUTPUT_PATH   "${GFP_SOURCE_DIR}")

# Set to build shared objects as option
OPTION(BUILD_SHARED_LIBS "Build package with shared libraries." ON)

# ---------------------------- Checking the OS ------------------------------
IF("${CMAKE_SYSTEM}" MATCHES "Linux")
	MESSAGE(STATUS "Compiling on Linux")
	SET(LINUX 1)
	EXEC_PROGRAM("uname" ARGS -m OUTPUT_VARIABLE myArch)
	SET(DIRARCHNAME 32)
	IF("${myArch}" MATCHES "x86_64")
		MESSAGE(STATUS "  64-bit architecture detected")
		SET(DIRARCHNAME 64)
	ENDIF("${myArch}" MATCHES "x86_64")
	SET(LINUX64 1)
	SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib${DIRARCHNAME}")
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib${DIRARCHNAME}")
	
ELSEIF(APPLE)
	SET(LINUX 0)
	MESSAGE(STATUS "Compiling on MacOSX")
ELSEIF(WIN32)
	SET(LINUX 0)
	IF(CYGWIN)
		MESSAGE(STATUS "Compiling on CygWin")
		SET(CYGWIN 1)
	ELSE(CYGWIN)
		MESSAGE(FATAL_ERROR "Compiling on Windows is not possible at the moment")
	ENDIF(CYGWIN)
ENDIF("${CMAKE_SYSTEM}" MATCHES "Linux")


# Library: Boost
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost 1.4.0 COMPONENTS program_options)
#FIND_PACKAGE(Boost COMPONENTS filesystem program_options thread unit_test_framework serialization)
IF(Boost_FOUND)
	MESSAGE(STATUS "Boost found")
	IF($ENV{VERBOSE})
 		MESSAGE(STATUS "   Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
 		MESSAGE(STATUS "   Boost_LINK_DIRECTORIES = ${Boost_LINK_DIRECTORIES}")
		MESSAGE(STATUS "   Boost Version = m:${Boost_MINOR_VERSION} M:${Boost_MAJOR_VERSION} v:${Boost_VERSION}")
	ENDIF($ENV{VERBOSE})
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
	LINK_DIRECTORIES(${Boost_LINK_DIRECTORIES})
	#ADD_DEFINITIONS(-DBOOST_VERSION=${Boost_VERSION})
ELSE(Boost_FOUND)
	MESSAGE(FATAL_ERROR "	Boost not found, it is REQUIRED to build the FLIRT.\n Boost can be found on http://www.boost.org")
ENDIF(Boost_FOUND)
 
 
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

# build

#Optimization
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -O3 -funroll-loops -fomit-frame-pointer -Wall -W -pipe -ftree-vectorize")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3 -funroll-loops -fomit-frame-pointer -Wall -W -pipe -ftree-vectorize")

#Debug (it can clash with optimization)
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3 -ggdb3 -fno-inline")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -ggdb3 -fno-inline")
ADD_DEFINITIONS(-Wall)


#Recursive compilation

ADD_SUBDIRECTORY(src)
