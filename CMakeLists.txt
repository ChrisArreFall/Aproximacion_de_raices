cmake_minimum_required(VERSION 3.10.2)


project(tarea03 VERSION 0.0.0 LANGUAGES CXX)


add_subdirectory (src)

enable_testing ()
add_subdirectory (test)
add_subdirectory(benchmarks)
