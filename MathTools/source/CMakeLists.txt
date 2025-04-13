set(SRC ${SRC}
${CMAKE_CURRENT_SOURCE_DIR}/ADT.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cellVolume.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cellVolume.f90
${CMAKE_CURRENT_SOURCE_DIR}/AxisSeperating.cpp
${CMAKE_CURRENT_SOURCE_DIR}/FileManager.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kaiser.cpp
${CMAKE_CURRENT_SOURCE_DIR}/LocalGlobalMap.cpp
${CMAKE_CURRENT_SOURCE_DIR}/math.cpp
${CMAKE_CURRENT_SOURCE_DIR}/median.cpp
${CMAKE_CURRENT_SOURCE_DIR}/OBBox.cpp
${CMAKE_CURRENT_SOURCE_DIR}/triangle.cpp
${CMAKE_CURRENT_SOURCE_DIR}/vector3d.cpp
)

set(HEADER ${HEADER}
${CMAKE_CURRENT_SOURCE_DIR}/ADT.h
${CMAKE_CURRENT_SOURCE_DIR}/AxisSeperating.h
${CMAKE_CURRENT_SOURCE_DIR}/cellNd.h
${CMAKE_CURRENT_SOURCE_DIR}/edgeNd.h
${CMAKE_CURRENT_SOURCE_DIR}/FileManager.h
${CMAKE_CURRENT_SOURCE_DIR}/LocalGlobalMap.h
${CMAKE_CURRENT_SOURCE_DIR}/OBBox.h
${CMAKE_CURRENT_SOURCE_DIR}/Operations.h
${CMAKE_CURRENT_SOURCE_DIR}/pointNd.h
${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
${CMAKE_CURRENT_SOURCE_DIR}/triangle.h
${CMAKE_CURRENT_SOURCE_DIR}/vector3d.h
${CMAKE_CURRENT_SOURCE_DIR}/vectorNd.h
)

add_library(MathTool STATIC ${SRC})

target_link_libraries(MathTool PUBLIC glm)



set(INCLUDE_DIRS include)
set(LIBRARIES MathTool)
set(LIB_DIR lib)

write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION 1.1.1
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION lib/cmake
    PATH_VARS INCLUDE_DIRS LIBRARIES LIB_DIR
    INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}

)

install(TARGETS MathTool DESTINATION lib)
install(FILES ${HEADER} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake)

if(BuildAsSubmodule)
    target_include_directories(MathTool PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()