cmake_minimum_required(VERSION 3.0.0)

add_executable(rar-engine nrm.cpp)
add_executable(bernoulli bernoulli.cpp)
target_include_directories(rar-engine PRIVATE "${CMAKE_ROOT}/../../include")
target_include_directories(bernoulli PRIVATE "${CMAKE_ROOT}/../../include")

target_link_libraries(bernoulli "${CMAKE_ROOT}/../../lib/libgsl.a")

# Compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17")

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -ggdb") # debug compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Werror")
endif()
