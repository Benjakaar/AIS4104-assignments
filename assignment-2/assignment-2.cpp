#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double deg)
{
    return deg * M_PI / 180;
}

double rad_to_deg(double rad)
{
    return rad * 180 / M_PI;
}


void example(double constant)
{
    Eigen::Matrix3d identity;
    identity <<
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    std::cout << "I: " << std::endl << identity << std::endl << std::endl;
    std::cout << constant <<"*I: " << std::endl << constant * identity << std::endl << std::endl;
}


Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& v)
{
    Eigen::Matrix3d SkewSymmetric;
    SkewSymmetric <<
        0.0, -v.z(), v.y(),
        v.z(), 0.0, -v.x(),
        -v.y(),v.x(), 0.0 ;
    return SkewSymmetric;
}

void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});

    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout <<-skew_matrix.transpose() << std::endl;
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                const Eigen::Vector3d &y,
                                                const Eigen::Vector3d &z) {
    Eigen::Matrix3d rotation_matrix;
    rotation_matrix <<
        x(0), y(0), z(0),
        x(1), y(1), z(1),
        x(2), y(2), z(2);

    return rotation_matrix;
}

Eigen::Matrix3d rotate_x(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        1.0, 0.0, 0.0,
        0.0, std::cos(radians), -std::sin(radians),
        0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_y(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        std::cos(radians), 0.0, std::sin(radians),
        0.0, 1, 0.0,
        -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_z(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        std::cos(radians), -std::sin(radians), 0.0,
        std::sin(radians), std::cos(radians), 0.0,
        0.0,0.0, 1.0;
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees) {
    double theta = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        std::cos(theta) + axis.x() * axis.x() * (1-std::cos(theta)), axis.y() * axis.z() * (1-std::cos(theta)) - axis.z()*std::sin(theta), axis.x() * axis.z() * (1-std::cos(theta)) - axis.y()*std::sin(theta),
        axis.x() * axis.y() * (1-std::cos(theta)) + axis.z()*std::sin(theta), std::cos(theta) + axis.y() * axis.y() * (1-std::cos(theta)), axis.y() * axis.z() * (1-std::cos(theta)) - axis.x()*std::sin(theta),
        axis.x() * axis.z() * (1-std::cos(theta)) - axis.y()*std::sin(theta), axis.y() * axis.z() * (1-std::cos(theta)) + axis.x()*std::sin(theta), std::cos(theta) + axis.z() * axis.z() * (1-std::cos(theta));
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    Eigen::Matrix3d euler_matrix;
    double x = deg_to_rad(e.x());
    double y = deg_to_rad(e.y());
    double z = deg_to_rad(e.z());

    double  cos_alpha = std::cos(x);
    double  cos_beta = std::cos(y);
    double  cos_gamma = std::cos(z);
    double  sin_alpha = std::sin(x);
    double  sin_beta = std::sin(y);
    double  sin_gamma = std::sin(z);

    euler_matrix <<
        cos_alpha * cos_beta, cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma, cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma,
        sin_alpha * cos_beta, sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma, sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma,
        -sin_beta, cos_beta * sin_gamma, cos_beta * cos_gamma;
    return euler_matrix;
}

void rotation_matrix_test()
{
    Eigen::Matrix3d rot = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0,-45.0, 90.0});
    Eigen::Matrix3d rot_aa = rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa = rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107}, Eigen::Vector3d{-0.5,-0.5, 0.707107}, Eigen::Vector3d{0.707107,-0.707107, 0.0});
    std::cout << "Rotation matrix from Euler: " << std::endl; std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl; std::cout << "Rotation matrix from frame axes: " << std::endl; std::cout << rot_fa << std::endl << std::endl;
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &d) {
    Eigen::Matrix4d matrix;
    matrix <<
        r(0,0), r(0,1), r(0,2), d(0),
        r(1,0), r(1,1), r(1,2), d(1),
        r(2,0), r(2,1), r(2,2), d(2),
        0.0,0.0,0.0,1.0;
    return matrix;
}

void transformation_matrix_test()
{
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45,-45.0, 90.0});
    Eigen::Vector3d v{1.0,-2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl; std::cout << transformation_matrix(r, v) << std::endl;


}

void transform_vector() {
    double euler_x = 60.0;  //X-rotation
    double euler_y = 45.0;  //Y-rotation
    double euler_z = 0;     //Z-rotation

    //Translation of {a} along the z-axis
    Eigen::Vector3d p(0.0, 0.0, 10.0);

    //Vector v_a
    Eigen::Vector3d v_a (2.5, 3.0, -10.0);

    //Rotation matrix from the Euler angels
    Eigen::Matrix3d rotation_matrix = rotation_matrix_from_euler_zyx(Eigen::Vector3d{euler_x,euler_y, euler_z});

    // Transformation matrix using the rotation and translation
    Eigen::Matrix4d trans_mat = transformation_matrix(rotation_matrix, p);

    //v_a as a 4x1 homogenous vector
    Eigen::Vector4d v_a_homogeneous (v_a.x(), v_a.y(), v_a.z(), 1);

    // Transforming v_a to v_w
    Eigen::Vector4d v_w_homogeneous = trans_mat * v_a_homogeneous;

    //Extrating the x, y, and z coodinates
    Eigen::Vector3d v_w (v_w_homogeneous.x(), v_w_homogeneous.y(), v_w_homogeneous.z());

    // Printing the answer
    std::cout << "Transformed vector in frame {w}: ["<< v_w.x()<<","<< v_w.y()<<"," <<v_w.z() <<"]"<< std::endl;
}

int main()
{
    //skew_symmetric_test();
    //rotation_matrix_test();
    //transformation_matrix_test();
    transform_vector();
    return 0;
}