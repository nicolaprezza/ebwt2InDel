cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (ebwt2snp)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(filter_snp filter_snp.cpp)
add_executable(vcf_vs_vcf vcf_vs_vcf.cpp)
add_executable(sam2vcf sam2vcf.cpp)
add_executable(snp2fastq snp2fastq.cpp)
add_executable(clust2snp clust2snp.cpp)
add_executable(ebwt2clust ebwt2clust.cpp)
add_executable(snp_vs_vcf snp_vs_vcf.cpp)
add_executable(differentialVCF differentialVCF.cpp)
