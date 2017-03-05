#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  for(unsigned int i =0; i< estimations.size(); i++){
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    float px2_plus_py2 = (px * px) + (py * py);
    float sqrt_px2_plus_py2 = sqrt(px2_plus_py2);
    float sqrt_px2_plus_py2_pow3 = sqrt_px2_plus_py2 * sqrt_px2_plus_py2 * sqrt_px2_plus_py2;
    //check division by zero
    if(px2_plus_py2 > 0.0001  && sqrt_px2_plus_py2 > 0.0001 && sqrt_px2_plus_py2_pow3 > 0.0001){
        float Hj_0_0 = (px/sqrt_px2_plus_py2);
        float Hj_0_1 = (py / sqrt_px2_plus_py2);
        Hj << Hj_0_0, Hj_0_1, 0, 0,
              (-1 * py / px2_plus_py2),(px / px2_plus_py2), 0, 0,
              ((py *((vx*py) - (vy*px)))/sqrt_px2_plus_py2_pow3 ), ((px *((vy*px) - (vx*py)))/sqrt_px2_plus_py2_pow3), Hj_0_0, Hj_0_1;
    }
    else cout<< "ERROR!! divide by zero\n";
    return Hj;
}