#Version definition
cmake_minimum_required(VERSION 3.20)
project(LungCancerIdentification)

#C++ standard definition
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
list(APPEND CMAKE_PREFIX_PATH "/usr/local/Halide/")

add_executable(LungCancerIdentification main.cpp)

#Include Halide library
find_package(Halide REQUIRED)
target_link_libraries(LungCancerIdentification PRIVATE Halide::Halide)

#Include dicom reader library
set(DICOM_USE_DCMTK ON)
