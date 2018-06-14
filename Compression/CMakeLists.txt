
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
 vtkMultiresolutionIOSeq.cxx
 vtkMultiresolutionIOSeqII.cxx
 vtkSparseMatrix.cxx
 vtkWaveletSubdivisionFilter.cxx
 )
 

# Set the added libraries
SET(LIB_ADDED vtkRangeCoding vtkSurface ${VTK_LIBRARIES})

TARGET_LINK_LIBRARIES(vtkWavemesh ${LIB_ADDED})

SET(WAVEMESH_EXAMPLES
SWA_TWA
wavemesh
DisplayWavelet
)
FOREACH(loop_var ${WAVEMESH_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkWavemesh ${LIB_ADDED})
ENDFOREACH(loop_var)

