#include <iostream>
#include <Eigen/Dense>

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &d) {
    Eigen::Matrix4d matrix;
    matrix <<
        r(0,0), r(0,1), r(0,2), d(0),
        r(1,0), r(1,1), r(1,2), d(1),
        r(2,0), r(2,1), r(2,2), d(2),
        0.0,0.0,0.0,1.0;
    return matrix;
}

double deg_to_rad(double deg)
{
    return deg * M_PI / 180;
}

double rad_to_deg(double rad)
{
    return rad * 180 / M_PI;
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

Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r) { //SIDE 498 I BOKA VISER FORMELENE
    Eigen::Vector3d angles;
    double beta, alpha, gamma;

    if (abs(r(2,0)) != 1) {
        beta = atan2(-r(2,0), sqrt(r(0,0)*r(0,0) + r(1,0)*r(1,0)));
        alpha = atan2(r(1,0), r(0,0));
        gamma = atan2(r(2,1), r(2,2));
    }

    else if (r(2,0) == -1){
        beta = std::numbers::pi/2;
        alpha = 0;
        gamma = atan2(r(0,1), r(1,1));
    }
    else{
        beta = -std::numbers::pi/2;
        alpha = 0;
        gamma = -atan2(r(0,1), r(1,1));
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

//--------------------------- Task 2b) ---------------------------
void Example3_28() {
    Eigen::VectorXd F_h(6);
    F_h << 0.0, 0.0, 0.0, 0.0, -5.0, 0.0;

    Eigen::VectorXd F_a(6);
    F_a << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix4d T_hf;
    T_hf <<
        1.0, 0.0, 0.0, -0.1,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix4d T_af;
    T_af <<
        1.0, 0.0, 0.0, -0.25,
        0.0, 0.0, 1.0, 0.0,
        0.0, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    Eigen::MatrixXd AD_Thf_transpose(6,6);
    AD_Thf_transpose  = adjoint_matrix(T_hf).transpose();

    Eigen::MatrixXd AD_Taf_transpose(6,6);
    AD_Taf_transpose  = adjoint_matrix(T_af).transpose();

    Eigen::VectorXd F_f(6);
    F_f = AD_Thf_transpose * F_h + AD_Taf_transpose * F_a;
    std::cout << F_f << std::endl;
}

//--------------------------- Task 3a) ---------------------------
Eigen::Matrix3d matrix_exponential_rot(const Eigen::Vector3d &w, double theta) {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d W = skew_symmetric(w);
    Eigen::Matrix3d R = I + std::sin(theta) * W + (1 - std::cos(theta)) * W * W;
    return R;
}

//--------------------------- Task 3b) ---------------------------
std::pair<Eigen::Vector3d, double> matrix_logarithm_rot(const Eigen::Matrix3d &r) {
    Eigen::Matrix3d skew_w;
    double theta;

    Eigen::Vector3d w;
    double r_trace = r.trace();

    if (r.isApprox(Eigen::Matrix3d::Identity())) {
        theta = 0;
        w << 0.0, 0.0, 0.0;
    }
    else if (r_trace == -1) {
        theta = M_PI;
        Eigen::Vector3d multi = {r(1,3), r(2,3), 1 + r(3,3)};
        w = (1/std::sqrt(2*(1+r(3,3))))*multi;
    }
    else {
        theta = std::acos(0.5*(r_trace-1));
        skew_w = (1/(2*std::sin(theta)))*(r-r.transpose());
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
    Eigen::Matrix3d R;
    Eigen::Vector3d p;
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double theta;
    Eigen::VectorXd exponential_coords(6);


    R = t.block(0, 0, 3, 3);
    p = t.block(0, 3, 3, 1);

    if (R.isApprox(I)) {
        w = {0,0,0};
        v = p/p.norm();
        theta = deg_to_rad(p.norm());
    }
    else {
        std::make_pair(w, theta) = matrix_logarithm_rot(R);
        theta = deg_to_rad(theta);
    }
    exponential_coords << v, w;
    return std::make_pair(exponential_coords, theta);
}


//--------------------------- Task 4a) ---------------------------
void print_pose(const std::string &label, const Eigen::Matrix4d &tf) {
    Eigen::Matrix3d R = tf.block(0, 0, 3, 3);
    Eigen::Vector3d p = tf.block(0, 3, 3, 1);

    Eigen::Vector3d euler_ZYX = euler_zyx_from_rotation_matrix(R);

    std::cout << "Position (x, y, z): ["
              << p.x() << ", "
              << p.y() << ", "
              << p.z() << "]" << std::endl;

    std::cout << "Orientation (rZ, rY, rX) in radians: ["
              << euler_ZYX[0] << ", "
              << euler_ZYX[1] << ", "
              << euler_ZYX[2] << "]"
              << std::endl;
}

//--------------------------- Task 4b) ---------------------------
Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions) {
    const double L1 = 10.0;
    const double L2 = 10.0;
    const double L3 = 10.0;

    double theta_1 = deg_to_rad(joint_positions[0]);
    double theta_2 = deg_to_rad(joint_positions[1]);
    double theta_3 = deg_to_rad(joint_positions[2]);

    Eigen::Matrix4d T_01, T_12, T_23, T_34, T_04;
    T_01 <<
        std::cos(theta_1), -std::sin(theta_1),0,0,
        std::sin(theta_1), std::cos(theta_1), 0,0,
        0,0,1,0,
        0,0,0,1;


    T_12 <<
        std::cos(theta_2), -std::sin(theta_2),0,L1,
        std::sin(theta_2), std::cos(theta_2), 0,0,
        0,0,1,0,
        0,0,0,1;

    T_23 <<
        std::cos(theta_3), -std::sin(theta_3),0,L2,
        std::sin(theta_3), std::cos(theta_3), 0,0,
        0,0,1,0,
        0,0,0,1;

    T_34 <<
        1,0,0,L3,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;

    T_04 = T_01*T_12*T_23*T_34;

    return T_04;
}

//--------------------------- Task 4c) ---------------------------
Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions) {
    const double L1 = 10.0;
    const double L2 = 10.0;
    const double L3 = 10.0;

    double theta_1 = deg_to_rad(joint_positions[0]);
    double theta_2 = deg_to_rad(joint_positions[1]);
    double theta_3 = deg_to_rad(joint_positions[2]);

    Eigen::Matrix4d M, T_04, skew_S_1, skew_S_2, skew_S_3;

    M <<
        1, 0, 0, L1 + L2 + L3,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

    Eigen::Vector3d w1(0, 0, 1);
    Eigen::Vector3d v1;
    v1 <<  0,0,0;

    Eigen::Vector3d w2(0, 0, 1);
    Eigen::Vector3d v2;
    v2 <<  0, -L1, 0;

    Eigen::Vector3d w3(0, 0, 1);
    Eigen::Vector3d v3;
    v3 << 0, -(L1 + L2), 0;

    Eigen::Matrix4d T1 = matrix_exponential_trans(w1, v1, theta_1);
    Eigen::Matrix4d T2 = matrix_exponential_trans(w2, v2, theta_2);
    Eigen::Matrix4d T3 = matrix_exponential_trans(w3, v3, theta_3);

    T_04 = T1 * T2 * T3 * M;

    return T_04;
}

//--------------------------- Task 5a) ---------------------------
Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions) {
    Eigen::Matrix4d M, T_06;
    const double H1 = 0.152;   //d1
    const double L1 = 0.244;   //a2
    const double L2 = 0.213;   //a3
    const double W1 = 0.131;   //d4
    const double H2 = 0.085;   //d5
    const double W2 = 0.092;   //d6

    double theta_1 = deg_to_rad(joint_positions[0]);
    double theta_2 = deg_to_rad(joint_positions[1]);
    double theta_3 = deg_to_rad(joint_positions[2]);
    double theta_4 = deg_to_rad(joint_positions[3]);
    double theta_5 = deg_to_rad(joint_positions[4]);
    double theta_6 = deg_to_rad(joint_positions[5]);

    // q, s, h
    Eigen::Vector3d q1(0, 0, 0), s1(0, 0, 1);  // Joint 1 (rotates about z-axis)
    Eigen::Vector3d q2(0, 0, H1), s2(0, 1, 0); // Joint 2 (rotates about y-axis)
    Eigen::Vector3d q3(L1, 0, H1), s3(0, 1, 0); // Joint 3 (rotates about y-axis)
    Eigen::Vector3d q4(L1 + L2, 0, H1), s4(0, 1, 0); // Joint 4 (rotates about y-axis)
    Eigen::Vector3d q5(L1 + L2, W1, H1), s5(0, 0, 1); // Joint 5 (rotates about z-axis)
    Eigen::Vector3d q6(L1 + L2, W1 + W2, H1), s6(0, 1, 0); // Joint 6 (rotates about y-axis)

    double h = 0;  // pitch = 0

    Eigen::VectorXd S1 = screw_axis(q1, s1, h);
    Eigen::VectorXd S2 = screw_axis(q2, s2, h);
    Eigen::VectorXd S3 = screw_axis(q3, s3, h);
    Eigen::VectorXd S4 = screw_axis(q4, s4, h);
    Eigen::VectorXd S5 = screw_axis(q5, s5, h);
    Eigen::VectorXd S6 = screw_axis(q6, s6, h);

    M <<
        -1, 0, 0, L1 + L2,
        0, 0, 1, W1 + W2,
        0, 1, 0, H1 - H2,
        0, 0, 0, 1;

    Eigen::Matrix4d T1 = matrix_exponential_trans({S1(0), S1(1), S1(2)}, {S1(3), S1(4), S1(5)}, theta_1);
    Eigen::Matrix4d T2 = matrix_exponential_trans({S2(0), S2(1), S2(2)}, {S2(3), S2(4), S2(5)}, theta_2);
    Eigen::Matrix4d T3 = matrix_exponential_trans({S3(0), S3(1), S3(2)}, {S3(3), S3(4), S3(5)}, theta_3);
    Eigen::Matrix4d T4 = matrix_exponential_trans({S4(0), S4(1), S4(2)}, {S4(3), S4(4), S4(5)}, theta_4);
    Eigen::Matrix4d T5 = matrix_exponential_trans({S5(0), S5(1), S5(2)}, {S5(3), S5(4), S5(5)}, theta_5);
    Eigen::Matrix4d T6 = matrix_exponential_trans({S6(0), S6(1), S6(2)}, {S6(3), S6(4), S6(5)}, theta_6);

    T_06 = T1 * T2 * T3 * T4 * T5 * T6 * M;

    return T_06;
}

//--------------------------- Task 5b) ---------------------------

Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions) {
    Eigen::Matrix4d  T_01, T_12, T_23, T_34, T_45, T_56, T_06, T_6END;
    Eigen::Matrix3d R01, R12, R23, R34, R45, R56;
    Eigen::Matrix3d R6END;
    Eigen::Vector3d P01, P12, P23, P34, P45, P56, PEND;

    const double H1 = 0.152;   //d0
    const double L1 = 0.244;   //a1
    const double L2 = 0.213;   //a2
    const double W1 = 0.131;   //d3
    const double H2 = 0.085;   //d4
    const double W2 = 0.092;   //d5

    R01 = rotate_z(joint_positions[0]);
    R12 = rotate_y(joint_positions[1]);
    R23 = rotate_y(joint_positions[2]);
    R34 = rotate_y(joint_positions[3]);
    R45 = rotate_y(-180)*rotate_z(joint_positions[4]); // Add a rotation of -180 degrees around y to tranform the frame
    R56 = rotate_x(-90)*rotate_z(joint_positions[5]); // Add a rotation of -90 degrees to around x tranform the frame
    R6END << rotate_z(0); // Identity matrix since there is no additional rotation between joint 6 and the end effector

    P01 = {0, 0, 0};
    P12 = {0, 0, H1};
    P23 = {L1, 0, 0};
    P34 = {L2, 0, 0};
    P45 = {0, W1, 0};
    P56 = {0, W2, 0};
    PEND = {0, -H2, 0};

    T_01 = transformation_matrix(R01, P01);
    T_12 = transformation_matrix(R12, P12);
    T_23 = transformation_matrix(R23, P23);
    T_34 = transformation_matrix(R34, P34);
    T_45 = transformation_matrix(R45, P45);
    T_56 = transformation_matrix(R56, P56);
    T_6END = transformation_matrix(R6END, PEND);

    T_06 =  T_01 * T_12 *T_23 * T_34 * T_45 * T_56 * T_6END;

    return T_06;
}

int main()
{
    std::vector<std::vector<double>> test_joint_positions = {
        {0.0, 0.0, 0.0},        // j1
        {90, 0.0, 0.0},         // j2
        {0.0, 90, 0.0},         // j3
        {0.0, 0.0, 90},         // j4
        {10.0, -15.0, 2.75 }    // j5
    };

    Eigen::Matrix4d T1 = planar_3r_fk_transform(test_joint_positions[2]);
    Eigen::Matrix4d T2 = planar_3r_fk_screw(test_joint_positions[2]);
    print_pose("Pose Description:", T1);
    print_pose("Pose Description:", T2);

    std::vector<std::vector<double>> test_joint_positions_6d = {
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},           // j1
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},         // j2
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}             // j3
    };

    Eigen::Matrix4d ur3_T_S = ur3e_fk_screw(test_joint_positions_6d[2]);
    Eigen::Matrix4d ur3_T_T = ur3e_fk_transform(test_joint_positions_6d[2]);

    //print_pose("Pose Description:", ur3_T_S);
    //print_pose("Pose Description:", ur3_T_T);

    return 0;
}