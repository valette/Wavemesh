 include_directories(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
)
 
 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(vtkSurface
 RenderWindow.cxx
 vtkDelaunay.cxx
 vtkDualMeshDisplay.cxx
 vtkOFFReader.cxx
 vtkOFFWriter.cxx
 vtkQuadricTools.cxx
 vtkRandomTriangulation.cxx
 vtkSMFReader.cxx
 vtkSurface.cxx
 vtkSurfaceBase.cxx
 vtkVolumeProperties.cxx
 vtkNeighbourhoodComputation.cxx
 vtkCurvatureMeasure.cxx
 )

target_link_libraries(vtkSurface ${VTK_LIBRARIES})

SET(COMMON_EXAMPLES
viewer
viewer2
RandomTriangulation
ExamplevtkSurface
)
FOREACH(loop_var ${COMMON_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkSurface ${VTK_LIBRARIES})
ENDFOREACH(loop_var)
