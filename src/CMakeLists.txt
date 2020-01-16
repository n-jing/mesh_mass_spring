FIND_PACKAGE(Eigen3 REQUIRED)
LIST(APPEND CMAKE_MODULE_PATH /home/jing/Workspace/Project/libigl/cmake)

INCLUDE(libigl)

SET(src_mesh
  solver.cpp
  writer.cpp
  mesh.cpp
  remove_duplicate_vert.cpp)
ADD_LIBRARY(mesh ${src_mesh})
TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC ${CMAKE_SOURCE_DIR}/inc)
TARGET_LINK_LIBRARIES(mesh PUBLIC Eigen3::Eigen igl::core bfgs)

ADD_EXECUTABLE(mass_spring main.cpp)
TARGET_LINK_LIBRARIES(mass_spring mesh)

SET(EXECUTABLE_OUTPUT_PATH "../bin")
