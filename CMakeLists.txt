cmake_minimum_required(VERSION 3.20)
project(arap_interactive CXX)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# Eigen
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(eigen)

# libigl
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG        v2.5.0
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(libigl)

add_subdirectory(src)
add_subdirectory(tests)