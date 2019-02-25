#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "math.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;
/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Laser measurement matrix
  H_laser_  << 1, 0, 0, 0,
               0, 1, 0, 0;

  // state covariance matrix P
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

  // the initial transition matrix F_
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

// Modify the F matrix so that the time is integrated
void FusionEKF::updateF(float time_elapsed)
{
  ekf_.F_(0, 2) = time_elapsed;
  ekf_.F_(1, 3) = time_elapsed;
}

// Modify the Q matrix
void FusionEKF::updateQ(float time_elapsed)
{
  float dt   = time_elapsed;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  // set the acceleration noise components
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    cout << "Kalman Filter Initialization " << endl;
    ekf_.x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      ekf_.x_ = tools.ConvertPolarToCartesian(measurement_pack.raw_measurements_[0], //ro
                                              measurement_pack.raw_measurements_[1], //theta
                                              measurement_pack.raw_measurements_[2]); //ro_dot
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
                 measurement_pack.raw_measurements_[1], 
                 0, 
                 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  
  // Set the previous timestamp with the current timestamp
  previous_timestamp_ = measurement_pack.timestamp_;
  
  //Update the state transition matrix with elapsed time in seconds
  updateF(dt);
  //Update the process covariance matrix Q
  updateQ(dt);
  
  /**
   * Prediction
   */
  ekf_.Predict();
  
  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_  = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_  = R_radar_;  
    ekf_.Hx_ = tools.PredictRadarMeasurement(ekf_.x_[0],
                                             ekf_.x_[1],
                                             ekf_.x_[2],
                                             ekf_.x_[3]); 
  } 
  else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Hx_ = ekf_.H_* ekf_.x_;  
  }
  ekf_.Update(measurement_pack.raw_measurements_);
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
