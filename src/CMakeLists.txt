FIND_PACKAGE(Eigen3 REQUIRED)
LIST(APPEND CMAKE_MODULE_PATH /home/jing/Workspace/Project/libigl/cmake)

INCLUDE(libigl)

# ADD_DEFINITIONS(-DRANDOM) #random vert weight and random k
# ADD_DEFINITIONS(-DRIGID)  #rigid mesh

SET(src_mesh
  solver.cpp
  writer.cpp
  mesh.cpp
  time_integral.cpp
  remove_duplicate_vert.cpp)
ADD_LIBRARY(mesh ${src_mesh})
TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC ${CMAKE_SOURCE_DIR}/inc)
TARGET_LINK_LIBRARIES(mesh PUBLIC Eigen3::Eigen igl::core bfgs)

ADD_EXECUTABLE(balance balance.cpp)
TARGET_LINK_LIBRARIES(balance mesh)

ADD_EXECUTABLE(dynamic dynamic.cpp)
TARGET_LINK_LIBRARIES(dynamic mesh)

SET(EXECUTABLE_OUTPUT_PATH "../bin")
