<!-- Start commands -->

cd /globalpath/to/arap-interactive
mkdir -p build && cd build
cmake .. && cmake --build . && ./tests/arap_tests
