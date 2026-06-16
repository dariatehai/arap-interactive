# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.13.0
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(googletest)

# Test executable
add_executable(arap_tests
    test_mesh.cpp
    test_ply_parser.cpp
    test_sparse_builder.cpp
)

target_include_directories(arap_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(arap_tests
    arap_lib
    gtest
    gtest_main
    Eigen3::Eigen
)

add_test(NAME ArnapTests COMMAND arap_tests)
