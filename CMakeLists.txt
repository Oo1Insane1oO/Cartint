cmake_minimum_required(VERSION 3.15.0)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

set(TESTING ON)

set(SOURCE source)
set(API interface)
set(INCLUDE include)
set(TEST test)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -lpthread")

project(
    Cartint
    LANGUAGES CXX
    DESCRIPTION "Solver for integrals over Gauss-Hermite functions in cartesian coordinates"
)

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES false ${SOURCE}/*.h ${SOURCE}/*.cpp ${API}/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${INCLUDE}/cartint.hpp)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE ${Boost_INCLUDE_DIRS} ${SOURCE}
    PUBLIC ${INCLUDE} ${API}
)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES} ${MPI_C})

if(${TESTING})
    set(CMAKE_BUILD_TYPE Debug)

    file(GLOB_RECURSE SOURCES LIST_DIRECTORIES false ${SOURCE}/*.h, ${SOURCE}/*.cpp)

    include(CTest)
    enable_testing()

    list(
        APPEND
        TESTNAMES
    )
    foreach(TESTNAME IN LISTS TESTNAMES)
        set(CARTINTSOURCES ${SOURCES} ${TEST}/${TESTNAME}.cpp)

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
