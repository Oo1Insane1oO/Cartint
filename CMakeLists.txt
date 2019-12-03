cmake_minimum_required(VERSION 3.15.0)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(TESTING OFF)

set(SOURCE source)
set(API interface)
set(INCLUDE include)
set(TEST test)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -lpthread -fPIC")

project(
    Cartint
    LANGUAGES CXX
    DESCRIPTION "Solver for integrals over Gauss-Hermite functions in cartesian coordinates"
)

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES false ${SOURCE}/*.h ${SOURCE}/*.hpp ${SOURCE}/*.cpp ${API}/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${INCLUDE}/cartint.hpp binding/bindings.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
find_package(PythonLibs 3 REQUIRED)

message("${PYTHON_INCLUDE_DIRS}")
return()

execute_process(
    COMMAND poetry show -v
    COMMAND grep virtualenv
    COMMAND /bin/sh -c "awk '{split($0, a); print a[3]}'"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE PYBIND_INCLUDE_DIRS
)
STRING(STRIP ${PYBIND_INCLUDE_DIRS} PYBIND_INCLUDE_DIRS)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
        ${Boost_INCLUDE_DIRS}
        ${PYBIND_INCLUDE_DIRS}/include/site/python3.7
        ${SOURCE}
        ${PYTHON_INCLUDE_DIRS}
    PUBLIC ${INCLUDE} ${API}
)
target_link_libraries(
    ${PROJECT_NAME} Eigen3::Eigen
        ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES} ${MPI_C} ${PYTHON_LIBRARIES}
)

if(${TESTING})
    set(CMAKE_BUILD_TYPE Debug)

    file(GLOB_RECURSE SOURCES LIST_DIRECTORIES false ${SOURCE}/*.h, ${SOURCE}/*.hpp ${SOURCE}/*.cpp binding/bindings.cpp)

    include(CTest)
    enable_testing()

    list(
        APPEND
        TESTNAMES
    )
    foreach(TESTNAME IN LISTS TESTNAMES)
        set(CARTINTSOURCES ${SOURCES} ${TEST}/${TESTNAME}.cpp binding/bindings.cpp)

        add_executable(${TESTNAME} ${CARTINTSOURCES})
        target_compile_features(${TESTNAME} PRIVATE cxx_std_17)

        find_package(Catch2 REQUIRED)
        target_include_directories(
            ${TESTNAME}
            PRIVATE ${SOURCE} ${TEST} ${ARMADILLO_INCLUDE_DIRS} ${CATCH2_INCLUDE_DIRS}
            PUBLIC ${INCLUDE}
        )
        target_link_libraries(
            ${TESTNAME}
            Eigen3::Eigen ${CATCH2_LIBRARIES} ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
        )
        
        add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:${TESTNAME}>)
    endforeach()
endif()
