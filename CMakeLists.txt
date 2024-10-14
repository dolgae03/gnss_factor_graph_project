cmake_minimum_required(VERSION 3.5)

project(ceres_solver_gnss)

find_package(Ceres REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options filesystem)
find_package(absl REQUIRED)
find_package(Eigen3 REQUIRED)

enable_testing()
add_test(NAME run_tests COMMAND run_tests)

# Boost include 디렉토리 설정
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(ceres_solver_gnss 
  src/main.cpp
  include/csv_utils.cpp
  include/cord_utils.cpp)
target_link_libraries(ceres_solver_gnss Ceres::ceres  ${Boost_LIBRARIES})