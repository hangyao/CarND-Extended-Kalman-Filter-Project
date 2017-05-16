#include <iostream>
#include <cmath>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cout << "Invalid input." << endl;
    return rmse;
  }

  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  MatrixXd Hj(3,4);
  Hj.fill(0.0);

  float px2_plus_py2 = px*px+py*py;
  if (px2_plus_py2 == 0) {
    px2_plus_py2 = 0.00000001;
  }

  float position_vec_magnitude = sqrt(px2_plus_py2);
  float position_vec_magnitude_cubed = pow(position_vec_magnitude, 3);

  float drho_dpx = px / position_vec_magnitude;
  float drho_dpy = py / position_vec_magnitude;
  float dphi_dpx = - py / px2_plus_py2;
  float dphi_dpy = px / px2_plus_py2;
  float drhodot_dpx = py * (vx * py - vy * px) / position_vec_magnitude_cubed;
  float drhodot_dpy = px * (vy * px - vx * py) / position_vec_magnitude_cubed;
  float drhodot_dvx = drho_dpx;
  float drhodot_dvy = drho_dpy;

  Hj << drho_dpx,    drho_dpy,    0,    0,
        dphi_dpx,    dphi_dpy,    0,    0,
        drhodot_dpx, drhodot_dpy, drhodot_dvx, drhodot_dvy;

  return Hj;
}

VectorXd Tools::PolarToCartesianMeasurement(const VectorXd& x_state) {
  float rho = x_state(0);
  float phi = x_state(1);
  float rho_dot = x_state(2);

  float tan_phi = tan(phi);
  float px = sqrt(rho*rho / (1 + tan_phi*tan_phi));
  float py = tan_phi * px;
  float vx = rho_dot * cos(phi);
  float vy = rho_dot * sin(phi);

  VectorXd cartesian_vec(4);
  cartesian_vec << px, py, vx, vy;
  return cartesian_vec;
}

VectorXd Tools::CartesianToPolarMeasurement(const VectorXd& x_state) {
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float position_vec_magnitude = sqrt(px*px+py*py);

  VectorXd polar_vec(3);
  polar_vec << 0,0,0;

  if (px == 0 || position_vec_magnitude == 0) {
    cout << "Division by 0 avoided." << endl;
    px = 0.00000001;
    position_vec_magnitude = sqrt(px*px+py*py);
  } else if (position_vec_magnitude == 0) {
    cout << "Division by 0 avoided." << endl;
    position_vec_magnitude = 0.00000001;
  }

  polar_vec << position_vec_magnitude,
               atan2(py, px),
               (px*vx+py*vy)/position_vec_magnitude;

  return polar_vec;
}
