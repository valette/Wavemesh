PROJECT(Wavemesh)

cmake_minimum_required(VERSION 3.5)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

OPTION(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# find_package
find_package(VTK
    REQUIRED
    COMPONENTS
    CommonSystem
    FiltersModeling
    ImagingGeneral
    InteractionStyle
    IOXML
    IOExport
    IOGeometry
    IOImport
    IOLegacy
    IOMINC
    IOPLY
    RenderingLabel
    RenderingOpenGL2
    RenderingGL2PSOpenGL2
)

include_directories(${VTK_INCLUDE_DIRS})

SET(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
SET(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/ACVD/Common)
SET(VTKWAVEMESH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Compression)
SET(VTKRANGECODING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/RangeCoding)

# The main library
ADD_SUBDIRECTORY(ACVD/Common)
ADD_SUBDIRECTORY(RangeCoding)
ADD_SUBDIRECTORY(Compression)
