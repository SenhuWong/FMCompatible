cmake_minimum_required(VERSION 3.16.0)

project(MathTool)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
# enable_language(C Fortran)
set(SRC)
set(HEADER)

include(CMakePackageConfigHelpers)
if(BuildAsSubmodule)

else()
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/glm)
endif()

add_subdirectory(source)

