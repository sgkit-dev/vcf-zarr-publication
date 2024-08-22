# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.1...3.29)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
  ts_afdist
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)


# Adding something we can run - Output name matches target name
add_executable(ts_afdist ts_afdist.cpp)

include(FetchContent)

FetchContent_Declare(
    tensorstore
    URL "https://github.com/google/tensorstore/archive/refs/tags/v0.1.64.tar.gz"
    URL_HASH SHA256=1dc632c6e9f83c033a2a16fe1a67ca38825902083ea0efe314769eee980c3126
)

# Additional FetchContent_Declare calls as needed...

FetchContent_MakeAvailable(tensorstore)

# Define a target that depends on TensorStore...

target_link_libraries(
    ts_afdist
    PRIVATE
    tensorstore::tensorstore tensorstore::all_drivers
    )
