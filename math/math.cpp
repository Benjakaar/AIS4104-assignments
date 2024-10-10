#include <iostream>
#include <Eigen/Dense>
#include "include/math.h"

namespace math {

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &d) {
        Eigen::Matrix4d matrix;
        matrix <<
            r(0,0), r(0,1), r(0,2), d(0),
            r(1,0), r(1,1), r(1,2), d(1),
            r(2,0), r(2,1), r(2,2), d(2),
            0.0,0.0,0.0,1.0;
        return matrix;
    }


    double deg_to_rads(double deg)
    {
        return deg * deg_to_rad_const;
    }

    double rad_to_degs(double rad)
    {
        return rad * rad_to_deg_const;
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

    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r) {
        Eigen::Vector3d angles;
        double beta, alpha, gamma;

        if (abs(r(2,0)) != 1) {
            beta = std::atan2(-r(2,0), std::sqrt(r(0,0)*r(0,0) + r(1,0)*r(1,0)));
            alpha = std::atan2(r(1,0), r(0,0));
            gamma = std::atan2(r(2,1), r(2,2));
        }

        else if (r(2,0) == -1){
            beta = std::numbers::pi/2;
            alpha = 0;
            gamma = std::atan2(r(0,1), r(1,1));
        }
        else{
            beta = -std::numbers::pi/2;
            alpha = 0;
            gamma = -std::atan2(r(0,1), r(1,1));
        }

        angles(0) = alpha;  //Z
        angles(1) = beta;   //Y
        angles(2) = gamma;  //X

        return angles;
        }

    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v) {
        //[ws, vs]
        Eigen::VectorXd Twist;
        Twist <<
            w(0),w(1),w(2),v(0),v(1),v(2);
        return Twist;
    }

    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h) {
        Eigen::VectorXd Twist(6);
        double theta_dot = 1.0;

        Eigen::Vector3d w = s * theta_dot;

        Eigen::Vector3d v = (-s*theta_dot).cross(q) + h*s*theta_dot;

        Twist <<
            w(0),w(1),w(2),v(0),v(1),v(2);
        return Twist;
    }

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf) {
        Eigen::Matrix3d R;
        Eigen::Vector3d p;
        Eigen::MatrixXd Adj(6,6);

        R = tf.block(0, 0, 3, 3);
        p = tf.block(0, 3, 3, 1);

        Eigen::Matrix3d skew_p = skew_symmetric(p);
        Eigen::Matrix3d pR = skew_p * R;

        Adj <<
            R(0,0), R(0,1), R(0,2), 0, 0, 0,
            R(1,0), R(1,1), R(1,2), 0, 0, 0,
            R(2,0), R(2,1), R(2,2), 0, 0, 0,
            pR(0,0), pR(0,1), pR(0,2),R(0,0),R(0,1),R(0,2),
            pR(1,0), pR(1,1), pR(1,2),R(1,0),R(1,1),R(1,2),
            pR(2,0), pR(2,1), pR(2,2),R(2,0),R(2,1),R(2,2);


        return Adj;
    }

    double cot(double x) {
        return 1.0 / std::tan(x);
    }

    //--------------------------- Task 2a) ---------------------------
    Eigen::Matrix3d rotate_x(double radians)
    {
        Eigen::Matrix3d matrix;
        matrix <<
            1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
        return matrix;
    }

    Eigen::Matrix3d rotate_y(double radians)
    {
        Eigen::Matrix3d matrix;
        matrix <<
            std::cos(radians), 0.0, std::sin(radians),
            0.0, 1, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);
        return matrix;
    }

    Eigen::Matrix3d rotate_z(double radians)
    {
        Eigen::Matrix3d matrix;
        matrix <<
            std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0,0.0, 1.0;
        return matrix;
    }

    void wrench_w_and_s() {
        Eigen::Vector3d f_w;
        Eigen::Vector3d m_s;
        Eigen::Vector3d e_ws;

        f_w << -30.0, 0.0, 0.0;
        m_s << 0.0, 0.0, 2.0;
        e_ws << 60.0, -60.0, 0.0; //in Euler YZX

        Eigen::Matrix3d R_ws = rotate_y(e_ws(0))*rotate_z(e_ws(1))*rotate_x(e_ws(2));
        Eigen::Vector3d f_s = R_ws*f_w;
        Eigen::Vector3d m_w = R_ws*m_s;

        std::cout << "f_w: ["<< f_w.x()<<","<< f_w.y()<<"," <<f_w.z() <<"]"<< std::endl;
        std::cout << "m_w: ["<< m_w.x()<<","<< m_w.y()<<"," <<m_w.z() <<"]"<< std::endl;
        std::cout << "f_s: ["<< f_s.x()<<","<< f_s.y()<<"," <<f_s.z() <<"]"<<"  <---flipped?"<< std::endl;
        std::cout << "m_s: ["<< m_s.x()<<","<< m_s.y()<<"," <<m_s.z() <<"]"<< std::endl;

    }

    //--------------------------- Task 3a) ---------------------------
    Eigen::Matrix3d matrix_exponential_rot(const Eigen::Vector3d &w, double theta) {
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d W = skew_symmetric(w);
        Eigen::Matrix3d R = I + std::sin(theta) * W + (1.0 - std::cos(theta)) * W * W;
        return R;
    }

    //--------------------------- Task 3b) ---------------------------
    std::pair<Eigen::Vector3d, double> matrix_logarithm_rot(const Eigen::Matrix3d &r) {
        Eigen::Matrix3d skew_w;
        double theta;
        Eigen::Vector3d w;
        double r_trace = r.trace();

        if (r.isApprox(Eigen::Matrix3d::Identity())) {
            theta = 0.0;
            w << 0.0, 0.0, 0.0;
        }
        else if (r_trace == -1.0) {
            theta = M_PI;
            Eigen::Vector3d multi = {r(0,2), r(1,2), 1 + r(2,2)};
            w = (1/std::sqrt(2*(1.0+r(2,2))))*multi;
        }
        else {
            theta = std::acos(0.5*(r_trace-1.0));
            skew_w = (1.0/(2*std::sin(theta)))*(r-r.transpose());
            w << skew_w(1,2), skew_w(2,0), skew_w(0,1);
        }
        return std::make_pair(w, theta);
    }
    //--------------------------- Task 3c) ---------------------------
    Eigen::Matrix4d matrix_exponential_trans(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta) {
        Eigen::Matrix4d T;
        Eigen::Matrix3d skew_w = skew_symmetric(w);
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d R = matrix_exponential_rot(w, theta);

        Eigen::Vector3d star = (I*theta + (1 - std::cos(theta)) * skew_w + (theta - std::sin(theta)) * skew_w * skew_w)*v;

        if (w.norm()== 1){
        T <<
            R(0,0),R(0,1), R(0,2), star.x(),
            R(1,0), R(1,1), R(1,2), star.y(),
            R(2,0), R(2,1), R(2,2), star.z(),
            0,0,0,1;
        }
        else {
            T <<
                I(0,0), I(0,1), I(0,2), v(0)*theta,
                I(1,0), I(1,1), I(1,2), v(1)*theta,
                I(2,0), I(2,1), I(2,2), v(2)*theta,
                0, 0, 0, 1;
        }
        return T;

    }

    //--------------------------- Task 3d) ---------------------------
    std::pair<Eigen::VectorXd, double> matrix_logarithm_trans(const Eigen::Matrix4d &t) {
        // ||p|| is the vector norm, aka the length.
        Eigen::Matrix3d R = t.block(0, 0, 3, 3);
        Eigen::Vector3d p = t.block(0, 3, 3, 1);
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

        Eigen::Vector3d w;
        Eigen::Vector3d v;
        double theta;
        Eigen::VectorXd exponential_coords(6);

        if (R.isApprox(I)) {
            w = Eigen::Vector3d::Zero();
            v = p/p.norm();
            theta = p.norm();
        }
        else {
            std::tie(w, theta) = matrix_logarithm_rot(R);  // Get unit vector w and angle theta
            Eigen::Matrix3d w_hat = skew_symmetric(w);
            v = (1.0 / theta *I - 0.5 * w_hat + (1.0 / theta - 0.5 * cot(theta / 2.0)) * w_hat * w_hat) * p;
        }
        exponential_coords << w, v;
        return std::make_pair(exponential_coords, theta);
    }


    //--------------------------- Task 4a) ---------------------------
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf) {
        Eigen::Matrix3d R = tf.block(0, 0, 3, 3);
        Eigen::Vector3d p = tf.block(0, 3, 3, 1);

        Eigen::Vector3d euler_ZYX = euler_zyx_from_rotation_matrix(R);

        std::cout << label << std::endl;

        std::cout << "Position (x, y, z): ["
                  << p.x() << ", "
                  << p.y() << ", "
                  << p.z() << "]" << std::endl;

        std::cout << "Orientation (rZ, rY, rX) in radians: ["
                  << euler_ZYX[0] << ", "
                  << euler_ZYX[1] << ", "
                  << euler_ZYX[2] << "]"
                  << std::endl;

        std::cout << std::endl;
    }
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
        Eigen::Matrix3d euler_matrix;
        double x = deg_to_rads(e.x());
        double y = deg_to_rads(e.y());
        double z = deg_to_rads(e.z());

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

}