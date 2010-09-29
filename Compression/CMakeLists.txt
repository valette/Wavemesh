
INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
    ${VTKRANGECODING_INCLUDE_DIR}
)

 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(vtkWavemesh
 vtkArithmeticCoder.cxx
 vtkMultiresolutionIO.cxx
 vtkSparseMatrix.cxx
 vtkWaveletSubdivisionFilter.cxx
 )
 

# Set the added libraries
SET(LIB_ADDED 
vtkRangeCoding
vtkSurface
vtkCommon
vtkFiltering
vtkImaging
vtkIO
vtkRendering
vtkHybrid
)

TARGET_LINK_LIBRARIES(vtkWavemesh ${LIB_ADDED})

SET(WAVEMESH_EXAMPLES
wavemesh
DisplayWavelet
)
FOREACH(loop_var ${WAVEMESH_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkSurface vtkWavemesh vtkRangeCoding
	vtkCommon vtkFiltering vtkGraphics vtkImaging vtkIO vtkRendering vtkHybrid)
ENDFOREACH(loop_var)

