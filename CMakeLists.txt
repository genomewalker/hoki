cmake_minimum_required(VERSION 3.20)
project(lhi CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# prefer conda-env packages
list(PREPEND CMAKE_PREFIX_PATH
    "$ENV{CONDA_PREFIX}"
    "/maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo")

find_package(zstd REQUIRED)

add_executable(lhi src/main.cpp)

find_path(XXHASH_INC xxhash.h
    HINTS "$ENV{CONDA_PREFIX}" "/maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo"
    PATH_SUFFIXES include REQUIRED)
find_library(XXHASH_LIB xxhash
    HINTS "$ENV{CONDA_PREFIX}" "/maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo"
    PATH_SUFFIXES lib REQUIRED)

target_include_directories(lhi PRIVATE src ${XXHASH_INC})
target_link_libraries(lhi PRIVATE ${XXHASH_LIB})

# zstd: handle both config-mode and find-module results
if(TARGET zstd::libzstd_shared)
    target_link_libraries(lhi PRIVATE zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
    target_link_libraries(lhi PRIVATE zstd::libzstd_static)
else()
    find_library(ZSTD_LIB zstd REQUIRED)
    target_link_libraries(lhi PRIVATE ${ZSTD_LIB})
    find_path(ZSTD_INC zstd.h REQUIRED)
    target_include_directories(lhi PRIVATE ${ZSTD_INC})
endif()

target_compile_options(lhi PRIVATE -O3 -march=native -Wall -Wextra -Wpedantic)

install(TARGETS lhi DESTINATION bin)
