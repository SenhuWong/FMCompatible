cmake_minimum_required(VERSION 3.10)

project(test)

set(BuildAsSubmodule ON)

set(CMAKE_CXX_STANDARD 20)

set(HDF5_DIR ~/thirdParty/HDF5/cmake)

find_package(MPI REQUIRED)

find_package(HDF5 REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/glm)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/Eigen)

add_subdirectory(MathTools)

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
${CMAKE_CURRENT_SOURCE_DIR}/StencilHandler.cpp
${CMAKE_CURRENT_SOURCE_DIR}/CellBrick.cpp)

add_executable(myAPP ${SRC})

target_include_directories(myAPP PUBLIC ${MPI_CXX_INCLUDE_DIRS})

target_link_libraries(myAPP PUBLIC Eigen3::Eigen MathTool ${MPI_LIBRARIES} ${HDF5_LIBRARIES})