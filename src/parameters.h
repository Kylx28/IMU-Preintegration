#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <eigen3/Eigen/Dense>

constexpr double ACC_N = 0.2;
constexpr double ACC_W = 0.02;
constexpr double GYR_N = 0.0002;
constexpr double GYR_W = 2.0e-5;

constexpr double O_P = 0;
constexpr double O_R = 3;
constexpr double O_V = 6;
constexpr double O_BA = 9;
constexpr double O_BG = 12;

constexpr double O_AN = 0;
constexpr double O_GN = 3;
constexpr double O_AW = 6;
constexpr double O_GW = 9;

// Realistic gravity vector (approximate Earth's gravity)
const Eigen::Vector3d G(0.0, 0.0, 9.81); // Assuming Earth's gravity is pointing downwards in the z-axis


#endif // PARAMETERS_H
