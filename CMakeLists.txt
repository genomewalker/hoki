cmake_minimum_required(VERSION 3.20)
project(hoki CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(PREPEND CMAKE_PREFIX_PATH
    "$ENV{CONDA_PREFIX}"
    "/maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo")

find_package(zstd REQUIRED)

add_executable(hoki src/main.cpp)
target_include_directories(hoki PRIVATE src)

if(TARGET zstd::libzstd_shared)
    target_link_libraries(hoki PRIVATE zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
    target_link_libraries(hoki PRIVATE zstd::libzstd_static)
else()
    find_library(ZSTD_LIB zstd REQUIRED)
    target_link_libraries(hoki PRIVATE ${ZSTD_LIB})
    find_path(ZSTD_INC zstd.h REQUIRED)
    target_include_directories(hoki PRIVATE ${ZSTD_INC})
endif()

target_compile_options(hoki PRIVATE -O3 -march=native -Wall -Wextra -Wpedantic)

install(TARGETS hoki DESTINATION bin)
