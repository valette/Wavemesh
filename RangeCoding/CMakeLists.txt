
INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIR}
)

 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(vtkRangeCoding
 QuasiStaticModel.cxx
 RangeDecoder.cxx
 RangeEncoder.cxx
 vtkArithmeticCoderBase.cxx
 )
 
TARGET_LINK_LIBRARIES(vtkRangeCoding vtkCommon)

