//
// Created by benka on 02.10.2024.
//

#ifndef MATH_H
#define MATH_H

namespace math {

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &d);

    const double deg_to_rad_const = 0.0174532925;

    const double rad_to_deg_const = 57.29578;

    double deg_to_rads(double deg);


    double rad_to_degs(double rad);

    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);

    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& v);


    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r);


    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);

    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);

    double cot(double x);

    Eigen::Matrix3d rotate_x(double degrees);

    Eigen::Matrix3d rotate_y(double degrees);

    Eigen::Matrix3d rotate_z(double degrees);

    Eigen::Matrix3d matrix_exponential_rot(const Eigen::Vector3d &w, double theta);

    std::pair<Eigen::Vector3d, double> matrix_logarithm_rot(const Eigen::Matrix3d &r);

    Eigen::Matrix4d matrix_exponential_trans(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);

    std::pair<Eigen::VectorXd, double> matrix_logarithm_trans(const Eigen::Matrix4d &t);

    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);

}
#endif //MATH_H
