cmake_minimum_required(VERSION 3.10)
project(SuccinctBitvector VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_TESTS "Build test programs" ON)
option(BUILD_BENCHMARKS "Build benchmark programs" ON)
option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(SUCCINCT_BUILD_BOTH "Build both static and shared libraries" ON)
option(SUCCINCT_HEADER_ONLY "Use header-only mode" OFF)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    if(WIN32 AND MINGW)
        # MinGW on Windows doesn't support sanitizers well
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    else()
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address,undefined")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /RTC1")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Check if source file exists
set(MAIN_SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/succinct_bitvector.cpp")
if(NOT EXISTS ${MAIN_SOURCE_FILE})
    message(FATAL_ERROR "Cannot find source file: ${MAIN_SOURCE_FILE}\n"
            "Please ensure the file exists at: src/succinct_bitvector.cpp\n"
            "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Create library
if(SUCCINCT_HEADER_ONLY)
    add_library(succinct_bitvector INTERFACE)
    target_include_directories(succinct_bitvector INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    )
    target_compile_features(succinct_bitvector INTERFACE cxx_std_17)
else()
    # Build both static and shared if requested
    if(SUCCINCT_BUILD_BOTH)
        # Static library
        add_library(succinct_bitvector_static STATIC src/succinct_bitvector.cpp)
        set_target_properties(succinct_bitvector_static PROPERTIES
                OUTPUT_NAME succinct_bitvector
                VERSION ${PROJECT_VERSION}
                POSITION_INDEPENDENT_CODE ON
        )
        target_include_directories(succinct_bitvector_static PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>
        )

        # Shared library
        add_library(succinct_bitvector_shared SHARED src/succinct_bitvector.cpp)
        set_target_properties(succinct_bitvector_shared PROPERTIES
                OUTPUT_NAME succinct_bitvector
                VERSION ${PROJECT_VERSION}
                SOVERSION 1
        )
        target_include_directories(succinct_bitvector_shared PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>
        )

        # Alias for default usage
        if(BUILD_SHARED_LIBS)
            add_library(succinct_bitvector ALIAS succinct_bitvector_shared)
        else()
            add_library(succinct_bitvector ALIAS succinct_bitvector_static)
        endif()
    else()
        # Single library based on BUILD_SHARED_LIBS
        add_library(succinct_bitvector src/succinct_bitvector.cpp)
        target_include_directories(succinct_bitvector PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>
        )
        set_target_properties(succinct_bitvector PROPERTIES
                VERSION ${PROJECT_VERSION}
                SOVERSION 1
                PUBLIC_HEADER include/succinct_bitvector.hpp
        )
    endif()
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_executable(test_succinct_bitvector test/test_succinct_bitvector.cpp)
    target_link_libraries(test_succinct_bitvector succinct_bitvector)
    add_test(NAME SuccinctBitvectorTest COMMAND test_succinct_bitvector)
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    add_executable(bench_succinct_bitvector bench/bench_succinct_bitvector.cpp)
    target_link_libraries(bench_succinct_bitvector succinct_bitvector)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_executable(bvtool examples/bvtool.cpp)
    target_link_libraries(bvtool succinct_bitvector)
endif()

# Installation
if(SUCCINCT_HEADER_ONLY)
    install(TARGETS succinct_bitvector
            EXPORT SuccinctBitvectorTargets
    )
    install(FILES
            include/succinct_bitvector.hpp
            include/succinct_bitvector_header_only.hpp
            DESTINATION include/succinct
    )
else()
    if(SUCCINCT_BUILD_BOTH)
        install(TARGETS succinct_bitvector_static succinct_bitvector_shared
                EXPORT SuccinctBitvectorTargets
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
                RUNTIME DESTINATION bin
        )
    else()
        install(TARGETS succinct_bitvector
                EXPORT SuccinctBitvectorTargets
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
                RUNTIME DESTINATION bin
                PUBLIC_HEADER DESTINATION include/succinct
        )
    endif()
    install(FILES include/succinct_bitvector.hpp
            DESTINATION include/succinct
    )
endif()

# Export targets
install(EXPORT SuccinctBitvectorTargets
        FILE SuccinctBitvectorTargets.cmake
        NAMESPACE Succinct::
        DESTINATION lib/cmake/SuccinctBitvector
)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/SuccinctBitvectorConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

# Create a simple config file if template doesn't exist
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SuccinctBitvectorConfig.cmake.in")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/SuccinctBitvectorConfig.cmake"
            "include(\"\${CMAKE_CURRENT_LIST_DIR}/SuccinctBitvectorTargets.cmake\")\n"
    )
else()
    configure_file(cmake/SuccinctBitvectorConfig.cmake.in
            "${CMAKE_CURRENT_BINARY_DIR}/SuccinctBitvectorConfig.cmake"
            @ONLY
    )
endif()

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SuccinctBitvectorConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SuccinctBitvectorConfigVersion.cmake"
        DESTINATION lib/cmake/SuccinctBitvector
)

# Documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in")
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        add_custom_target(docs
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM
        )
    endif()
endif()

# Summary
message(STATUS "")
message(STATUS "Configuration summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build tests: ${BUILD_TESTS}")
message(STATUS "  Build benchmarks: ${BUILD_BENCHMARKS}")
message(STATUS "  Build examples: ${BUILD_EXAMPLES}")
message(STATUS "  Build shared libs: ${BUILD_SHARED_LIBS}")
message(STATUS "")