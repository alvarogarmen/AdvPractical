cmake_minimum_required(VERSION 3.16)
project(AdvPractical)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mconsole")

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

set(SOURCES
        Edge.hpp
        Grader.hpp
        Graph.hpp
        InputGraphManually.hpp
        Node.hpp
        )
add_executable(advpractical ${SOURCES} main.cpp argtable3.c)

set_target_properties(advpractical PROPERTIES LINK_FLAGS "-mconsole")