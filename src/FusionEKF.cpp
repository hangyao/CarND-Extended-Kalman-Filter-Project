#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0,    0,
        0, 1, 0,    0,
        0, 0, 1000, 0,
        0, 0, 0,    1000;


  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x_(4);
    x_.fill(1.0);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ = tools.PolarToCartesianMeasurement(measurement_pack.raw_measurements_);
      Hj_ = tools.CalculateJacobian(x_);
      ekf_.Init(x_, P_, F_, Hj_, R_radar_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1],
            0, 0;
      ekf_.Init(x_, P_, F_, H_laser_, R_laser_);
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

   if (dt > 0.001) {
     // Modify the F matrix so that the time is integrated
     ekf_.F_(0,2) = dt;
     ekf_.F_(1,3) = dt;

     // Set the process covariance matrix Q
     float dt_2 = dt * dt;
     float dt_3 = dt_2 * dt;
     float dt_4 = dt_3 * dt;

     ekf_.Q_ = MatrixXd(4, 4);
     ekf_.Q_ << (dt_4*noise_ax/4), 0,                 (dt_3*noise_ax/2), 0,
                0,                 (dt_4*noise_ay/4), 0,                 (dt_3*noise_ay/2),
                (dt_3*noise_ax/2), 0,                 (dt_2*noise_ax),   0,
                0,                 (dt_3*noise_ay/2), 0,                 (dt_2*noise_ay);

    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  ekf_.z_pred_ = tools.CartesianToPolarMeasurement(ekf_.x_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.z_pred_ = ekf_.H_ * ekf_.x_;
  }

  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
