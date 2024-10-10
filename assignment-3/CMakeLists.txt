find_package(Eigen3 CONFIG REQUIRED)

add_executable(assignment-3 assignment-3.cpp)

target_link_libraries(assignment-3 PRIVATE Eigen3::Eigen)
target_link_libraries(assignment-3 PRIVATE math)