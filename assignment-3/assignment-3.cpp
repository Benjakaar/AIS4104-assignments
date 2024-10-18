#include <iostream>
#include <numeric>
#include <Eigen/Dense>
#include "../math/include/math.h"
#include <functional>


//  Task 1a)
Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v) {

    Eigen::VectorXd r(v.size());

    for(int i = 0; i < v.size(); i++) {
        r[i] = v[i];
    }

    return r;
}

// Task 1b)
bool is_average_below_eps(const std::vector<double> &values, double eps = 10e-7, uint8_t n_values = 5u) {

    if(values.size() < n_values) {
        return false;
    }

    const double sum = std::accumulate(values.end()-n_values, values.end(), 0.0);
    const bool check = std::abs(sum/n_values) < eps;

    return check;
}

// Task 1c)
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_space_chain() {
    double L1 = 0.2435;
    double L2 = 0.2132;
    double W1 = 0.1315;
    double W2 = 0.0921;
    double H1 = 0.1518;
    double H2 = 0.08535;

    Eigen::Matrix3d mr = math::rotate_y(-90.0 * math::deg_to_rad_const)
                        * math::rotate_x(-90.0 * math::deg_to_rad_const)
                        * math::rotate_z(-90.0 * math::deg_to_rad_const);

    Eigen::Matrix4d m = math::transformation_matrix(mr, Eigen::Vector3d{L1 + L2, W1 + W2, H1 - H2});


    std::vector<Eigen::VectorXd> screws{
        math::screw_axis({0, 0, 0}, {0, 0, 1}, 0),
        math::screw_axis({0, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1 + L2, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1 + L2, W1, 0}, {0, 0, -1}, 0),
        math::screw_axis({L1 + L2, 0, H1 - H2}, {0, 1, 0}, 0)
    };

    return std::make_pair(m, screws);
}

// Task 1d)
Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta) {
    return math::matrix_exponential_trans(screw.head<3>(), screw.tail<3>(), theta);
}


Eigen::Matrix4d ur3_space_fk(const Eigen::VectorXd &joint_positions) {
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();

    for(int i = 0; i < joint_positions.size(); i++) {

        t06 *= matrix_exponential(space_screws[i], joint_positions[i]);
    }
    return t06*m;
}

// Task 1e)
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain() {
        double L1 = 0.2435;
        double L2 = 0.2132;
        double W1 = 0.1315;
        double W2 = 0.0921;
        double H1 = 0.1518;
        double H2 = 0.08535;

    Eigen::Matrix3d mr = math::rotate_y(-90.0 * math::deg_to_rad_const)
                        * math::rotate_x(-90.0 * math::deg_to_rad_const)
                        * math::rotate_z(-90.0 * math::deg_to_rad_const);

    Eigen::Matrix4d m = math::transformation_matrix(mr, Eigen::Vector3d{L1 + L2, W1 + W2, H1 - H2});


    std::vector<Eigen::VectorXd> screws{
        math::screw_axis({0, 0, 0}, {0, 0, 1}, 0),
        math::screw_axis({0, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1 + L2, 0, H1}, {0, 1, 0}, 0),
        math::screw_axis({L1 + L2, W1, 0}, {0, 0, -1}, 0),
        math::screw_axis({L1 + L2, 0, H1 - H2}, {0, 1, 0}, 0)
    };


    std::vector<Eigen::VectorXd> body{
        math::adjoint_matrix(m.inverse())*screws[0],
        math::adjoint_matrix(m.inverse())*screws[1],
        math::adjoint_matrix(m.inverse())*screws[2],
        math::adjoint_matrix(m.inverse())*screws[3],
        math::adjoint_matrix(m.inverse())*screws[4],
        math::adjoint_matrix(m.inverse())*screws[5]
        };


    return std::make_pair(m, body);
}


Eigen::Matrix4d ur3_body_fk(const Eigen::VectorXd &joint_positions) {
    auto [m, body_screws] = ur3e_body_chain();
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();

    for(int i = 0; i < joint_positions.size(); i++) {

        t06 *= matrix_exponential(body_screws[i], joint_positions[i]);
    }
    return m*t06;
}

// Testing task 1
void ur3e_test_fk()
{
    std::cout << "Forward kinematics tests" << std::endl;

    math::print_pose("space transformation one: ",ur3_space_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));
    math::print_pose("body transformation one: ",ur3_body_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));
    std::cout << std::endl;
    math::print_pose("space transformation one: ",ur3_space_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0})*math::deg_to_rad_const));
    math::print_pose("body transformation one: ",ur3_body_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0})*math::deg_to_rad_const));
    std::cout << std::endl;
    math::print_pose("space transformation one: ",ur3_space_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));
    math::print_pose("body transformation one: ",ur3_body_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));
    std::cout << std::endl;
    math::print_pose("space transformation one: ",ur3_space_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));
    math::print_pose("body transformation one: ",ur3_body_fk(std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0})*math::deg_to_rad_const));

}


// Task 2a
std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double x_0, double dx_0 = 0.5, double eps = 10e-7) {
        int iterations = 0;
        double x_n = x_0;
        double f_x, df_x;

        while (true) {
            f_x = f(x_n);

            if (std::abs(f_x) < eps) {
                break;
            }

            df_x = (f(x_n + dx_0) - f(x_n)) / dx_0;  // Estimation of the derivative

            x_n = x_n - f_x / df_x;

            iterations++;
        }
        return std::make_pair(iterations, x_n);
    }

// Task 2b
std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f, double x_0, double gamma = 0.1, double dx_0 = 0.5, double eps = 10e-7) {
        int iterations = 0;
        double fn_1 = 0.0;
        double fn = f(x_0);
        double df = fn - fn_1;
        double x = x_0;
        double dx = dx_0;

        while (iterations < 1000) {
            if (std::abs(fn) < eps) { // If function is converged enough
                break;
            }
            fn_1 = fn;
            x += dx;
            fn = f(x);
            df = fn - fn_1;
            dx = -gamma * (dx/df) * fn;
            iterations++;

        }
        return std::make_pair(iterations, x);
    }


void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = newton_raphson_root_find(f, x0);
    std::cout << "NR root f, x0=" << x0 << "-> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) << std::endl;
}

void test_gradient_descent_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = gradient_descent_root_find(f, x0);
    std::cout << "GD root f, x0=" << x0 << "-> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) << std::endl;
}

void test_root_find() {
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x) { return (x- 3.0) * (x- 3.0)- 1.0; };
    test_newton_raphson_root_find(f1,-20);
    test_gradient_descent_root_find(f1,-20);
}

// Task 3a
Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions) {
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Identity(space_screws.size(), current_joint_positions.size());

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity(); // Store tranformation between loops

    for (int i = 0; i < current_joint_positions.size(); i++) {
        if (i == 0) {
            Jacobian.col(0) = space_screws[i]; // First index is just S1
        }
        else {
            T = T * (matrix_exponential(space_screws[i-1], current_joint_positions[i-1]));

            Eigen::MatrixXd Adjoint_T = math::adjoint_matrix(T);

            Jacobian.col(i) = Adjoint_T * space_screws[i];
        }
    }

    return Jacobian;
}

// Task 3b
Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions) {
    auto [m, body_screws] = ur3e_body_chain();
    Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Identity(body_screws.size(), current_joint_positions.size());

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity(); // Store tranformation between loops

    for (int i = current_joint_positions.size()-1; i >= 0; i--) {

        if (i == current_joint_positions.size() - 1) {
            Jacobian.col(i) = body_screws[i]; // Last index is just Bn
        }
        else {
            T = T * (matrix_exponential(-body_screws[i+1], current_joint_positions[i+1]));

            Eigen::MatrixXd Adjoint_T = math::adjoint_matrix(T);

            Jacobian.col(i) = Adjoint_T * body_screws[i];
        }
    }

    return Jacobian;
}

void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions)
{
    Eigen::Matrix4d tsb = ur3_body_fk(joint_positions);
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::MatrixXd jb = ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = math::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = math::adjoint_matrix(tsb.inverse());

    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js << std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb << std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb- ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js- ad_tsb * jb << std::endl << std::endl;

}

void ur3e_test_jacobian() {
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad_const);
    ur3e_test_jacobian(std_vector_to_eigen(std::vector<double>{45.0,-20.0, 10.0, 2.5, 30.0,-50.0}) * math::deg_to_rad_const);
}

// Task 4a
std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd
    &current_joint_positions, double gamma = 1e-3, double v_e = 4e-3, double w_e = 4e-3) {

    long iterations = 0;

    Eigen::VectorXd joint_positions = current_joint_positions;

    double norm_v;
    double norm_w;

    Eigen::VectorXd V_b(6);
    V_b << 0.0,0.0,0.0,0.0,0.0,0.0;
    double theta;

    while (iterations < 100000) {
        Eigen::Matrix4d t_sb = ur3_space_fk(joint_positions);

        std::tie(V_b,theta) = math::matrix_logarithm_trans(t_sb.inverse() * t_sd);

        Eigen::Vector3d w_b = V_b.head(3)*theta;
        Eigen::Vector3d v_b = V_b.tail(3)*theta;

        norm_v = v_b.norm();
        norm_w = w_b.norm();

        Eigen::MatrixXd J_b = ur3e_body_jacobian(joint_positions);
        Eigen::MatrixXd J_b_pinv = J_b.completeOrthogonalDecomposition().pseudoInverse();

        joint_positions += (-gamma * J_b_pinv * V_b);

        if (norm_v <= v_e && norm_w <= w_e) {
            break;
        }
        iterations++;
    }

    return std::make_pair(iterations, joint_positions);
}

//Task 4b
void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0) {
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = math::transformation_matrix(math::rotation_matrix_from_euler_zyx(zyx), pos);
    std::cout << t_sd << std::endl;
    auto [iterations, j_ik] = ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = ur3_body_fk(j_ik);
    math::print_pose(" IK pose",t_ik );
    math::print_pose("Desired pose",t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg_const << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg_const << std::endl << std::endl;
}

void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0) {
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = ur3_space_fk(joint_positions);
    auto [iterations, j_ik] = ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = ur3_body_fk(j_ik);
    math::print_pose(" IK pose",t_ik); math::print_pose("Desired pose",t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg_const << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * math::rad_to_deg_const << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg_const << std::endl << std::endl;
}

void ur3e_ik_test() {
    Eigen::VectorXd j_t0 = std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad_const;
    Eigen::VectorXd j_t1 = std_vector_to_eigen(std::vector<double>{0.0, 0.0,-89.0, 0.0, 0.0, 0.0}) * math::deg_to_rad_const;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0,-90.0} * math::deg_to_rad_const, j_t0);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0,-90.0} * math::deg_to_rad_const, j_t1);
    Eigen::VectorXd j_t2 = std_vector_to_eigen(std::vector<double>{50.0,-30.0, 20, 0.0,-30.0, 50.0}) * math::deg_to_rad_const;
    Eigen::VectorXd j_d1 = std_vector_to_eigen(std::vector<double>{45.0,-20.0, 10.0, 2.5, 30.0,-50.0}) * math::deg_to_rad_const;
    ur3e_ik_test_configuration(j_d1, j_t0); ur3e_ik_test_configuration(j_d1, j_t2);
}
int main() {
    ur3e_test_fk();
    test_root_find();
    ur3e_test_jacobian();
    ur3e_ik_test();
    return 0;
 }
