SET(hlbfgs
  HLBFGS_BLAS.cpp
  HLBFGS.cpp
  ICFS.cpp
  LineSearch.cpp
  Lite_Sparse_Matrix.cpp)

ADD_LIBRARY(bfgs ${hlbfgs})
TARGET_INCLUDE_DIRECTORIES(bfgs PUBLIC ${CMAKE_SOURCE_DIR}/thirdparty/HLBFGS)

