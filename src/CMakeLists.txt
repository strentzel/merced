cmake_minimum_required(VERSION 2.8.6)
add_executable(test convex.c cut.c path.c prism.c raster.c vec.c)
target_link_libraries(test m)

