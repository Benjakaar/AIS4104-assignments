find_package(Eigen3 CONFIG REQUIRED)

add_executable(assignment_1 assignment-1.cpp)

target_link_libraries(assignment_1 PRIVATE Eigen3::Eigen)

