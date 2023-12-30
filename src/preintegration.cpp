#include <eigen3/Eigen/Dense>
#include <vector>
#include "parameters.cpp"

/*
Preintegration algorithm:
Inputs - IMU acc and gyr measurements, dt, acceleration bias, gyroscope bias
1. Propagation from frame k to k+1
2. Update jacobian matrix
3. Store values in vector
4. If bias changes too much then repropagate, otherwise use first order approximation
*/

class Preintegration{
    private:
        //imu measurements, etc
        double dt, sum_dt;
        Eigen::Vector3d acc0, gyr0, acc1, gyr1, bias_a, bias_g;
        const Eigen::Vector3d inital_acc, inital_gyr;

        //Position, Velocity, and Orientation estiamtes
        Eigen::Vector3d delta_p, delta_v;
        Eigen::Quaterniond delta_q;

        //buffer vectors
        std::vector<double> dt_vec;
        std::vector<Eigen::Vector3d> acc_vec;
        std::vector<Eigen::Vector3d> gyr_vec;

        //Jacobian and Covariance matrices
        Eigen::Matrix<double, 15, 15> jacobian, covariance;

        //Noise Matrix
        Eigen::Matrix<double, 18, 18> noise;

    public:
        Preintegration() = delete;
        Preintegration(const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0, const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
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
            noise = Eigen::Matrix<double, 18, 18>::Zero();
            noise.block<3,3>(0,0) = (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
            noise.block<3,3>(3,3) = (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
            noise.block<3,3>(6,6) = (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
            noise.block<3,3>(9,9) = (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
            noise.block<3,3>(12,12) = (ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
            noise.block<3,3>(15, 15) = (GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
            
        }

        void eulerIntegration(double _dt, const Eigen::Vector3d &acc_0, const Eigen::Vector3d &gyr_0, const Eigen::Vector3d &_delta_p, const Eigen::Vector3d &_delta_v,
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

                //Configure Matrix F (I + F * dt)
                Eigen::Matrix<double,15,15> F;
                F = Eigen::Matrix<double,15,15>::Zero();
                //Row 1:
                F.block<3,3>(3,0) = Eigen::Matrix3d::Identity();
                //Row 2:
                F.block<3,3>(6,3) = -1 * rotationMatrix * acc_skew;
                F.block<3,3>(9,3) = -1 * rotationMatrix;
                //Row 3:
                F.block<3,3>(6,6) = -1 * gyr_skew;
                F.block<3,3>(12,6) = -1 * Eigen::Matrix3d::Identity();


                F = (Eigen::Matrix<double, 15, 15>::Identity() + (F * _dt));

                //Configure Matrix G (G * dt)
                Eigen::Matrix<double, 12, 15> G;
                G = Eigen::Matrix<double, 15, 12>::Zero();
                //Row 2 
                G.block<3,3>(3,0) = -1 * rotationMatrix;
                //Row 3
                G.block<3,3>(6,3) = -1 * Eigen::Matrix3d::Identity();
                //Row 4
                G.block<3,3>(9,6) = Eigen::Matrix3d::Identity();
                //Row 5
                G.block<3,3>(12,9) = Eigen::Matrix3d::Identity();

                //Covariance Matrix
                covariance = (F * covariance * F.transpose()) + (G * noise * G.transpose());

                //Jacobian Matrix
                jacobian = F * jacobian; //inital jacobian is the identity matrix
            }
        }

        void propagate(double _dt, const Eigen::Vector3d &acc_1, const Eigen::Vector3d &gyr_1){
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
        }
        
        void store_values(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr){
            dt_vec.push_back(dt);
            acc_vec.push_back(acc);
            gyr_vec.push_back(gyr);
            propagate(dt, acc, gyr);
        }

        void repropagate(const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
            sum_dt = 0;
            acc0 = inital_acc;
            gyr0 = inital_gyr;
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

        Eigen::Matrix<double, 15, 1>calc_residuals(const Eigen::Vector3d &Pi, const Eigen::Vector3d &Qi, const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                                  const Eigen::Vector3d &Pj, const Eigen::Vector3d &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj)
        {
            Eigen::Matrix<double, 15, 1> residuals;
            //Pre-integration term estimation
        }
};