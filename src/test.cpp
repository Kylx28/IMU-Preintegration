#include "preintegration.h"

using namespace Eigen;

int main(){
    Vector3d acc_0(1, 0, 0);
    Vector3d gyr_0(0, 1, 0);
    Vector3d bias_acc(0, 0, 0);
    Vector3d bias_gyr(0, 0, 0);

    Vector3d delta_p(0, 0, 0);
    Vector3d delta_v(0, 0, 0);
    Quaterniond delta_q(1, 0, 0, 0);

    double dt = 0.1;

    Preintegration p(acc_0, gyr_0, bias_acc, bias_gyr);

    p.propagate(dt, acc_0, gyr_0);
    
    return 0;
}
