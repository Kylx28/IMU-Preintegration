#include <eigen3/Eigen/Dense>

//these variables represent acceleration noise, acceleration bias, gyroscope noise, and gyroscope bias
//Normally taken from IMU config (pre-defined values)
double ACC_N = 0.2;
double ACC_W = 0.02;
double GYR_N = 0.0002;
double GYR_W = 2.0e-5;

//State Order and Noise Order (??)
double O_P = 0;
double O_R = 3;
double O_V = 6;
double O_BA = 9;
double O_BG = 12;

double O_AN = 0;
double O_GN = 3;
double O_AW = 6;
double O_GW = 9;

//Gravity Vector
Eigen::Vector3d G;

//void readParameters();
