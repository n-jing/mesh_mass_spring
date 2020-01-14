FIND_PACKAGE(Eigen3 REQUIRED)
LIST(APPEND CMAKE_MODULE_PATH /home/jing/Workspace/Project/libigl/cmake)

INCLUDE(libigl)

SET(src_mesh
  mesh.cpp)
ADD_LIBRARY(mesh ${src_mesh})
TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC ${CMAKE_SOURCE_DIR}/inc)
TARGET_LINK_LIBRARIES(mesh Eigen3::Eigen)


ADD_EXECUTABLE(mass_spring main.cpp)
TARGET_INCLUDE_DIRECTORIES(mass_spring PUBLIC ${CMAKE_SOURCE_DIR}/inc)
TARGET_LINK_LIBRARIES(mass_spring mesh igl::core)

SET(EXECUTABLE_OUTPUT_PATH "../bin")
