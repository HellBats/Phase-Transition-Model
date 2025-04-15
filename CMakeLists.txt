cmake_minimum_required(VERSION 3.28.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)


project(Perovskites VERSION 1.0)
find_package(OpenMP REQUIRED)
file(GLOB_RECURSE SRC_FILES src/*.c)
add_executable(build ${SRC_FILES})
target_link_libraries(build PRIVATE -lgsl)
target_link_libraries(build PRIVATE -lm)
add_compile_definitions(OMP_STACKSIZE=512M)


if(OpenMP_C_FOUND)
    target_link_libraries(build PRIVATE OpenMP::OpenMP_C)
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()
target_include_directories(build PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)