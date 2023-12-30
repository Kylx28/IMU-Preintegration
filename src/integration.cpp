#include <eigen3/Eigen/Dense>
#include <vector>

class Integration{
    private:
        double dt;
        double sum_dt;
        Eigen::Vector3d acc_0, gyr_0, acc_1, gyr_1, bias_a, bias_g;
        const Eigen::Vector3d _acc0, _gyr0; //store constant inital values
        //Covariance matrix declaration
        Eigen::Vector3d delta_p, delta_v; //propagated position and velocity constraints
        Eigen::Quaterniond delta_q;

        //buffer variables
        std::vector<double> dt_buf;
        std::vector<Eigen::Vector3d> acc_buf;
        std::vector<Eigen::Vector3d> gyr_buf;


    public:
        Integration() = delete;
        Integration(const Eigen::Vector3d &acc0, const Eigen::Vector3d &gyr0, const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
            acc_0 = acc0;
            gyr_0 = gyr0;
            bias_a = bias_acc;
            bias_g = bias_gyr;
            sum_dt = 0;
            //Initalize covariance
            //From VINS-Mono, at beginning:
            //position and velocity are 0, orientation is identity
            delta_p = Eigen::Vector3d::Zero();
            delta_q = Eigen::Quaterniond::Identity();
            delta_v = Eigen::Vector3d::Zero();

            //Implement noise matrix
        }

        void eulerIntegration(double _dt, const Eigen::Vector3d &acc0, const Eigen::Vector3d &gyr0, const Eigen::Vector3d &acc1, const Eigen::Vector3d &gyr1,
                              const Eigen::Vector3d &_delta_p, const Eigen::Vector3d &_delta_v, const Eigen::Quaterniond &_delta_q,
                              const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr,
                              Eigen::Vector3d &delta_p_1, Eigen::Vector3d &delta_v_1, Eigen::Quaterniond &delta_q_1,
                              Eigen::Vector3d &bias_acc_1, Eigen::Vector3d &bias_gyr_1, bool update_jacobian)
        {
            //Equation (7) terms
            delta_p_1 = _delta_p + (_delta_v * _dt) + (0.5 * (_delta_q * (acc0 - bias_acc)) * _dt * _dt); //alpha i+1
            delta_v_1 = _delta_v + (_delta_q * (acc0 - bias_acc) * _dt); //beta i + 1
            delta_q_1 = _delta_q * Eigen::Quaterniond(1, 0.5 * gyr0(0) * _dt, 0.5 * gyr0(1) * _dt, 0.5 * gyr0(2) * _dt); //gamma i + 1

            bias_acc_1 = bias_acc;
            bias_gyr_1 = bias_gyr;

            //code to update jacobian
            if(update_jacobian){
                //do something
            }
        }

        void propagate(double _dt, const Eigen::Vector3d &acc1, const Eigen::Vector3d &gyr1, int method){
            //propagate new estimates for integration factors (equation (7) in paper)
            dt = _dt;
            acc_1 = acc1;
            gyr_1 = gyr1;
            Eigen::Vector3d delta_p_1, delta_v_1, bias_acc_1, bias_gyr_1;
            Eigen::Quaterniond delta_q_1;

            if(method == 1){ //euler integration
                eulerIntegration(_dt, acc_0, gyr_0, acc_1, gyr_1, delta_p, delta_v, delta_q, bias_a, bias_g, delta_p_1, delta_v_1,
                                 delta_q_1, bias_acc_1, bias_gyr_1, 1);
            }
            else if(method == 2) { //mid point
                //do something
            }
            else if(method == 3){ //rk4
                //do something
            }

            delta_p = delta_p_1;
            delta_v = delta_v_1;
            delta_q = delta_q_1;
            bias_a = bias_acc_1;
            bias_g = bias_gyr_1;
            delta_q.normalize();
            sum_dt += dt;
            acc_0 = acc_1;
            gyr_0 = gyr_1;
        }

        void buffer(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr){
            dt_buf.push_back(dt);
            acc_buf.push_back(acc);
            gyr_buf.push_back(gyr);
            propagate(dt, acc, gyr, 1); //propagate with euler method
        }

        void repropagate(const Eigen::Vector3d &bias_acc, const Eigen::Vector3d &bias_gyr){
            //reset?
            sum_dt = 0;
            acc_0 = _acc0;
            gyr_0 = _gyr0;
            delta_p.setZero();
            delta_v.setZero();
            delta_q.setIdentity();
            bias_a = bias_acc;
            bias_g = bias_gyr;
            //set jacobian as well
            //set covariance
            //set buffer values in loop
        }
};