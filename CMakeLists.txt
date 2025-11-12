cmake_minimum_required(VERSION 3.10)
project(ProbabilityCalculator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable optimization flags for low-latency performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops")

# Add executable
add_executable(prob_calc main.cpp)

# Include header files
target_include_directories(prob_calc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Optional: Add warnings
target_compile_options(prob_calc PRIVATE -Wall -Wextra -Wpedantic)

# Link math library
target_link_libraries(prob_calc m)
