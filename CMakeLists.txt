PROJECT(Wavemesh)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
#MARK_AS_ADVANCED( FORCE CMAKE_BACKWARDS_COMPATIBILITY )

# for CMake 2.6 corrected behaviour (see "cmake --help-policy CMP0003")
IF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 NEW)
ENDIF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)
# --------------------------------------------------------------------------


# --------------------- OPTIONS OF THE PROJECT -----------------------------
# Build static lib by default

OPTION(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)
OPTION(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)

OPTION(BUILD_RANGECODING "Build Range Coder/Decoder." ON)
OPTION(BUILD_WAVEMESH "Build Wavemesh and Compression Tools." ON)
# Python is not working at this time, and not supported anymore unless someone requests it
#OPTION(BUILD_PYTHON "Build the Wrapped Python classes." OFF)
#
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# -------------------- DEPENDENCIES TO EXTERNALS ---------------------------
# Find and load VTK settings.

FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)


SET(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
SET(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Common)
SET(VTKWAVEMESH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Compression)
SET(VTKRANGECODING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/RangeCoding)
# --------------------------------------------------------------------------
# Find and load Doxygen settings (when required)


IF(BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
ENDIF(BUILD_DOCUMENTATION)

# --------------------------------------------------------------------
# Switch to enable Multithreaded clustering


# The main library
SUBDIRS(Common)

IF (BUILD_RANGECODING)
    SUBDIRS(RangeCoding)
ENDIF (BUILD_RANGECODING)

SUBDIRS(VectorCompression)


IF(BUILD_WAVEMESH)
	IF(NOT BUILD_RANGECODING)
	MESSAGE(FATAL_ERROR "You need to set RANGECODING build to ON in order to build Wavemesh")
	ENDIF(NOT BUILD_RANGECODING)
    SUBDIRS(Compression)
ENDIF(BUILD_WAVEMESH)


## eof - CMakeLists.txt
