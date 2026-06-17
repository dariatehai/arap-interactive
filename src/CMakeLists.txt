add_library(arap_lib
    mesh.cpp
    ply_parser.cpp
    sparse_builder.cpp
)

target_include_directories(arap_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(arap_lib PUBLIC
    Eigen3::Eigen
    igl::core
)

if(MSVC)
    target_compile_options(arap_lib PRIVATE /W4)
else()
    target_compile_options(arap_lib PRIVATE -Wall -Wextra)
endif()

add_executable(arap_app
    main.cpp
)

target_include_directories(arap_app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(arap_app PRIVATE
    arap_lib
    igl::glfw
)