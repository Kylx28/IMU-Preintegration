#ifndef PREINTEGRATION_H
#define PREINTEGRATION_H

#include <eigen3/Eigen/Dense>
#include <vector>
#include "parameters.h" // If needed
#include <iostream>

class Preintegration {
private:
    double dt, sum_dt;
    Eigen::Vector3d acc0, gyr0, acc1, gyr1, bias_a, bias_g;
    const Eigen::Vector3d initial_acc, initial_gyr;
    Eigen::Vector3d delta_p, delta_v;
    Eigen::Quaterniond delta_q;
    std::vector<double> dt_vec;
    std::vector<Eigen::Vector3d> acc_vec;
    std::vector<Eigen::Vector3d> gyr_vec;
    Eigen::Matrix<double, 15, 15> jacobian, covariance;
    Eigen::Matrix<double, 12, 12> noise;

public:
    Preintegration(const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0,
                   const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr);

    void eulerIntegration(double _dt, const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0,
                          const Eigen::Vector3d &_delta_p, const Eigen::Vector3d &_delta_v,
                          const Eigen::Quaterniond &_delta_q, const Eigen::Vector3d &bias_acc,
                          const Eigen::Vector3d &bias_gyr, Eigen::Vector3d &delta_p_1,
                          Eigen::Vector3d &delta_v_1, Eigen::Quaterniond &delta_q_1,
                          Eigen::Vector3d &bias_a_1, Eigen::Vector3d &bias_g_1, bool update_jacobian);

    void propagate(double _dt, const Eigen::Vector3d &acc_1, const Eigen::Vector3d &gyr_1);

    void store_values(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr);

    void repropagate(const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr);

    Eigen::Matrix<double, 15, 1> calc_residuals(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi,
                                                const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai,
                                                const Eigen::Vector3d &Bgi, const Eigen::Vector3d &Pj,
                                                const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj,
                                                const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj);
};

#endif  // PREINTEGRATION_H
