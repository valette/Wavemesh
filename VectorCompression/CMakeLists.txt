
INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIR}
    ${VTKRANGECODING_INCLUDE_DIR}
)


SET(VECTOR_COMPRESSION_PROGS
encodeCoeffs
decodeCoeffs
)
FOREACH(loop_var ${VECTOR_COMPRESSION_PROGS})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkRangeCoding vtkHybrid)
ENDFOREACH(loop_var)

