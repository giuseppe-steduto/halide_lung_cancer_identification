#Version definition
cmake_minimum_required(VERSION 3.20)
project(LungCancerIdentification)

#C++ standard definition
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
list(APPEND CMAKE_PREFIX_PATH "/usr/local/Halide/")

add_executable(LungCancerIdentification main.cpp roi_definition.cpp)

#Include Halide library
find_package(Halide REQUIRED)
#target_link_libraries(LungCancerIdentification PRIVATE Halide::Halide)

# search for DCMTK library and header files
find_package(DCMTK REQUIRED)
#Include dicom reader library
set(DICOM_USE_DCMTK ON)

# specify DCMTK header include directories
include_directories(${DCMTK_INCLUDE_DIRS})
find_library(DINIFTI_DCMTK_LIB NAMES dcmimgle PATHS ${LIBDIR})
target_link_libraries(LungCancerIdentification
        PRIVATE Halide::Halide
        dcmimgle
        jpeg
        png
        dcmdata
        ofstd
        )
