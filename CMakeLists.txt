######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.10)
project(simple_dotplot VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS )
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

#set(CMAKE_CXX_FLAGS "-fsanitize=address -fexceptions")  # leak checking and address sanitization
set(CMAKE_CXX_FLAGS "-fexceptions")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/Bam.cpp
        src/IterativeSummaryStats.cpp
        src/Sam.cpp
        )

project(wambam)
add_library(wambam STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(wambam PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(wambam PROPERTIES DEFINE_SYMBOL "")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- sparsepp ---------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

ExternalProject_Add(project_sparsepp
        GIT_REPOSITORY https://github.com/greg7mdp/sparsepp.git
        TIMEOUT 5
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "" # Disable install step, is a header only lib!
        )

add_dependencies(wambam project_sparsepp)

# Specify include dir
ExternalProject_Get_Property(project_sparsepp SOURCE_DIR)
include_directories(${SOURCE_DIR}/sparsepp)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- htslib ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# --------------------

message(STATUS "Fetching htslib")

# https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2

find_package(ZLIB REQUIRED)

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'curl' library")
endif()

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
)

FetchContent_MakeAvailable(project_htslib)

set(HTSLIB_SOURCES
        ${project_htslib_SOURCE_DIR}/cram/cram_codecs.c
        ${project_htslib_SOURCE_DIR}/cram/cram_decode.c
        ${project_htslib_SOURCE_DIR}/cram/cram_encode.c
        ${project_htslib_SOURCE_DIR}/cram/cram_external.c
        ${project_htslib_SOURCE_DIR}/cram/cram_index.c
        ${project_htslib_SOURCE_DIR}/cram/cram_io.c
        ${project_htslib_SOURCE_DIR}/cram/cram_samtools.c
        ${project_htslib_SOURCE_DIR}/cram/cram_stats.c
        ${project_htslib_SOURCE_DIR}/cram/files.c
        ${project_htslib_SOURCE_DIR}/cram/mFILE.c
        ${project_htslib_SOURCE_DIR}/cram/open_trace_file.c
        ${project_htslib_SOURCE_DIR}/cram/pooled_alloc.c
        ${project_htslib_SOURCE_DIR}/cram/rANS_static.c
        ${project_htslib_SOURCE_DIR}/cram/sam_header.c
        ${project_htslib_SOURCE_DIR}/cram/string_alloc.c
        ${project_htslib_SOURCE_DIR}/os/rand.c
        ${project_htslib_SOURCE_DIR}/bcf_sr_sort.c
        ${project_htslib_SOURCE_DIR}/bgzf.c
        ${project_htslib_SOURCE_DIR}/bgzip.c
        ${project_htslib_SOURCE_DIR}/errmod.c
        ${project_htslib_SOURCE_DIR}/faidx.c
        ${project_htslib_SOURCE_DIR}/hfile.c
        ${project_htslib_SOURCE_DIR}/hfile_gcs.c
        ${project_htslib_SOURCE_DIR}/hfile_libcurl.c
        ${project_htslib_SOURCE_DIR}/hfile_net.c
        ${project_htslib_SOURCE_DIR}/hts.c
        ${project_htslib_SOURCE_DIR}/htsfile.c
        ${project_htslib_SOURCE_DIR}/kfunc.c
        ${project_htslib_SOURCE_DIR}/knetfile.c
        ${project_htslib_SOURCE_DIR}/kstring.c
        ${project_htslib_SOURCE_DIR}/md5.c
        ${project_htslib_SOURCE_DIR}/multipart.c
        ${project_htslib_SOURCE_DIR}/probaln.c
        ${project_htslib_SOURCE_DIR}/realn.c
        ${project_htslib_SOURCE_DIR}/regidx.c
        ${project_htslib_SOURCE_DIR}/sam.c
        ${project_htslib_SOURCE_DIR}/synced_bcf_reader.c
        ${project_htslib_SOURCE_DIR}/tabix.c
        ${project_htslib_SOURCE_DIR}/tbx.c
        ${project_htslib_SOURCE_DIR}/textutils.c
        ${project_htslib_SOURCE_DIR}/thread_pool.c
        ${project_htslib_SOURCE_DIR}/vcf.c
        ${project_htslib_SOURCE_DIR}/vcf_sweep.c
        ${project_htslib_SOURCE_DIR}/vcfutils.c
        )

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")


set(HTS_DIR _deps/project_htslib-src)

# HTSLIB configure
add_custom_target(
        BUILD_HTS
        ALL
        #        OUTPUT ${PROJECT_SOURCE_DIR}/external/htslib/config.h # //TODO: at build time ${htslib_SOURCE_DIR} is not actually populated!! is an empty string
        WORKING_DIRECTORY ${HTS_DIR}
        COMMAND pwd
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
        COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
)


# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
        ${CMAKE_SOURCE_DIR}/external/
)

#add_custom_target(HTSLIB_CONFIGURED DEPENDS "${CMAKE_BINARY_DIR}/_deps/project_htslib-src/config.h")

add_library(htslib ${HTSLIB_SOURCES})
add_dependencies(htslib BUILD_HTS)
add_dependencies(wambam htslib)

target_link_libraries(htslib pthread bz2 z curl)

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################


# -------- TESTS --------

set(TESTS
        test_htslib_bam_reader
        )


foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            wambam
            htslib
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        wam
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            wambam
            htslib
            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        wambam
        DESTINATION ${BINARY_INSTALL_DIR})
