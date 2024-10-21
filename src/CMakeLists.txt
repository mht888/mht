cmake_minimum_required(VERSION 3.17)
project(myalice)


include_directories(${CMAKE_CURRENT_SOURCE_DIR} ../emp-tool /usr/include  ../PackedArray) #gmp还没下载，要在linux系统下载,这里是gmp的位置

find_package(SEAL 3.6 REQUIRED)
#find_package(OpenMP REQUIRED)
find_package(OpenSSL REQUIRED)
#find_package(GMP REQUIRED)

find_package(pybind11 REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(she-seal she-seal.cpp)

# 添加 pybind11 模块
pybind11_add_module(myalice myalice.cpp)
add_library(protocol SHARED myalice.cpp)

target_link_libraries(she-seal SEAL::seal OpenSSL::SSL)
#target_link_libraries(she-seal SEAL::seal OpenMP::OpenMP_CXX)



