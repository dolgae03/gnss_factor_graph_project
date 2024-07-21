cmake_minimum_required(VERSION 3.5)

project(ceres_solver_gnss)

find_package(Ceres REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(ceres_solver_gnss 
  src/main.cpp
  include/csv_utils.cpp
  include/cord_utils.cpp)
target_link_libraries(ceres_solver_gnss Ceres::ceres)