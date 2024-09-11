find_package(Eigen3 CONFIG REQUIRED)

add_executable(assignment_2 assignment-2.cpp)

target_link_libraries(assignment_2 PRIVATE Eigen3::Eigen)
