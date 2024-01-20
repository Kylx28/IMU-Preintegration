#include "preintegration.h"

Preintegration::Preintegration(const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0, const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
    acc0 = acc_0;
    gyr0 = gyr_0;
    bias_a = bias_acc;
    bias_g = bias_gyr;
    sum_dt = 0;

    //Initialize covariance and jacobian matrix
    covariance = Eigen::Matrix<double, 15, 15>::Zero();
    jacobian = Eigen::Matrix<double, 15, 15>::Identity();

    //Constraint Variables - Position, Velocity, Orientation
    delta_p = Eigen::Vector3d::Zero();
    delta_v = Eigen::Vector3d::Zero();
    delta_q = Eigen::Quaterniond::Identity();

    //Noise initalization
    noise = Eigen::Matrix<double, 12, 12>::Zero();
    noise.block<3,3>(0,0) = (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
    noise.block<3,3>(3,3) = (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
    noise.block<3,3>(6,6) = (ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
    noise.block<3,3>(9, 9) = (GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
    
}

void Preintegration::eulerIntegration(double _dt, const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0, const Eigen::Vector3d &_delta_p, const Eigen::Vector3d &_delta_v,
                        const Eigen::Quaterniond &_delta_q, const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr,
                        Eigen::Vector3d &delta_p_1, Eigen::Vector3d &delta_v_1, Eigen::Quaterniond &delta_q_1,
                        Eigen::Vector3d &bias_a_1, Eigen::Vector3d &bias_g_1, bool update_jacobian)
{
    Eigen::Vector3d temp_gyr = gyr_0 - bias_gyr;
    delta_p_1 = _delta_p + (_delta_v * _dt) + (0.5 * (_delta_q * (acc_0 - bias_acc) * _dt * _dt));
    delta_v_1 = _delta_v + (_delta_q * (acc_0 - bias_acc) * _dt);
    delta_q_1 = _delta_q * Eigen::Quaterniond(1, 0.5 * temp_gyr(0) * _dt, 0.5 * temp_gyr(1) * _dt, 0.5 * temp_gyr(2) * _dt);
    bias_a_1 = bias_acc;
    bias_g_1 = bias_gyr;

    //jacobian update
    if(update_jacobian){
        Eigen::Matrix3d rotationMatrix = _delta_q.toRotationMatrix();
        Eigen::Vector3d temp_acc = acc_0 - bias_acc;
        Eigen::Vector3d temp_gyr = gyr_0 - bias_gyr;
        Eigen::Matrix3d acc_skew, gyr_skew;
        acc_skew<<0, -1 * temp_acc(2), temp_acc(1), 
                    temp_acc(2), 0, -1 * temp_acc(0),
                    -1 * temp_acc(1), temp_acc(0), 0;
        gyr_skew<<0, -1 * temp_gyr(2), temp_gyr(1),
                    temp_gyr(2), 0, -1 * temp_gyr(0),
                    -1 * temp_gyr(1), temp_gyr(0), 0;

        //Configure Matrix F
        Eigen::Matrix<double,15,15> F;
        F = Eigen::Matrix<double,15,15>::Zero();
        F.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
        F.block<3,3>(0,3) = -0.5 * rotationMatrix * acc_skew * _dt;
        F.block<3,3>(0,6) = Eigen::Matrix3d::Identity() * _dt;
        F.block<3,3>(0,9) = -0.5 * gyr_skew * _dt * _dt;
        F.block<3,3>(3,3) = Eigen::Matrix3d::Identity() - gyr_skew * _dt;
        F.block<3,3>(3,12) = -1 * Eigen::Matrix3d::Identity() * _dt;
        F.block<3,3>(6,3) = -1 * rotationMatrix * acc_skew * _dt;
        F.block<3,3>(6,6) = Eigen::Matrix3d::Identity();
        F.block<3,3>(6,9) = -1 * rotationMatrix * _dt;
        F.block<3,3>(9,9) = Eigen::Matrix3d::Identity();
        F.block<3,3>(12,12) = Eigen::Matrix3d::Identity();


        F = (Eigen::Matrix<double, 15, 15>::Identity() + (F * _dt));

        //Configure Matrix G
        Eigen::Matrix<double, 15, 12> G;
        G = Eigen::Matrix<double, 15, 12>::Zero();
        G.block<3,3>(0,0) = -0.5 * rotationMatrix * _dt * _dt;
        G.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
        G.block<3,3>(6,0) = rotationMatrix;
        G.block<3,3>(9,6) = Eigen::Matrix3d::Identity();
        G.block<3,3>(12,12) = Eigen::Matrix3d::Identity();

        //Covariance Matrix
        covariance = (F * covariance * F.transpose()) + (G * noise * G.transpose());

        //Jacobian Matrix
        jacobian = F * jacobian; //initial jacobian is the identity matrix
    }
}

void Preintegration::propagate(double _dt, const Eigen::Vector3d &acc_1, const Eigen::Vector3d &gyr_1){
    dt = _dt;
    acc1 = acc_1;
    gyr1 = gyr_1;
    Eigen::Vector3d delta_p_1, delta_v_1, bias_a_1, bias_g_1;
    Eigen::Quaterniond delta_q_1;

    //Euler Implementation for now
    eulerIntegration(dt, acc0, gyr0, delta_p, delta_v, delta_q, bias_a, bias_g, delta_p_1, delta_v_1, delta_q_1, bias_a_1, bias_g_1, 1);

    delta_p = delta_p_1;
    delta_v = delta_v_1;
    delta_q = delta_q_1;
    delta_q.normalize();
    bias_a = bias_a_1;
    bias_g = bias_g_1;
    sum_dt += dt;
    acc0 = acc1;
    gyr0 = gyr1;

    std::cout << "delta_p:" << std::endl << delta_p << std::endl;
    std::cout << "delta_v:" << std::endl << delta_v << std::endl;
    std::cout << "delta_q:" << std::endl << delta_q.w() << std::endl << delta_q.vec() << std::endl;

}

void Preintegration::store_values(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr){
    dt_vec.push_back(dt);
    acc_vec.push_back(acc);
    gyr_vec.push_back(gyr);
    propagate(dt, acc, gyr);
}

void Preintegration::repropagate(const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
    sum_dt = 0;
    acc0 = initial_acc;
    gyr0 = initial_gyr;
    delta_p.setZero();
    delta_v.setZero();
    delta_q.setIdentity();
    bias_a = bias_acc;
    bias_g = bias_gyr;

    //set jacobian
    //set covariance

    for(int i = 0; i < dt_vec.size(); i++){
        propagate(dt_vec[i], acc_vec[i], gyr_vec[i]);
    }
}

Eigen::Matrix<double, 15, 1> Preintegration::calc_residuals(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                            const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj)
{
    Eigen::Matrix<double, 15, 1> residuals;

    Eigen::Matrix3d J_ba_p = jacobian.block<3,3>(O_P, O_BA);
    Eigen::Matrix3d J_bw_p = jacobian.block<3,3>(O_P, O_BG);
    Eigen::Matrix3d J_ba_v = jacobian.block<3,3>(O_V, O_BA);
    Eigen::Matrix3d J_bw_v = jacobian.block<3,3>(O_V, O_BG);
    Eigen::Matrix3d J_bw_q = jacobian.block<3,3>(O_R, O_BG);

    Eigen::Vector3d d_ba = Bai - bias_a;
    Eigen::Vector3d d_bw = Bgi - bias_g;

    Eigen::Vector3d temp_vec = J_bw_q * d_bw;

    //Pre-integration term estimation
    Eigen::Vector3d _delta_p = delta_p + J_ba_p * d_ba + J_bw_p * d_bw;
    Eigen::Vector3d _delta_v = delta_v + J_ba_v * d_ba + J_bw_p * d_bw;
    Eigen::Quaterniond _delta_q = delta_q * Eigen::Quaterniond(1, 0.5 * temp_vec(0), 0.5 * temp_vec(1), 0.5 * temp_vec(2));

    residuals.block<3,1>(0, 0) = Qi.inverse() * (Pj - Pi + 0.5*(G * dt * dt) - (Vi * dt)) - _delta_p;
    residuals.block<3,1>(3, 0) = Qi.inverse() * (Vj + (G * dt) - Vi) - _delta_v;
    residuals.block<3,1>(6, 0) = (_delta_q.inverse() * Qi.inverse() * Qj).vec();  //vector part of quaternion
    residuals.block<3,1>(9, 0) = Baj - Bai;
    residuals.block<3,1>(12, 0) = Bgj - Bgi;

    return residuals;
}
