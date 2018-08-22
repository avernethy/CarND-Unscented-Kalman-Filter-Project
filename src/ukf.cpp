#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30; // <- MODIFY THIS

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30; // <- MODIFY THIS
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
 long long previous_time_stamp;
 if (!is_initialized_) {
   x_ << 0,0,0,0,0;
   P_ << 0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << 0,0,0,0,0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      x_ << 0,0,0,0,0;
    }

    //update the previous time stamp
    previous_time_stamp = meas_package.timestamp_;

    //done initializing, no need to predict or update
    is_initialized_ = true;
    return;
 }

 double delta_t = (meas_package.timestamp_ - previous_time_stamp)/1000000.0;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  cout << delta_t << endl;
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_x;

  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  //set first colume of sigma point matrix
  Xsig.col(0)  = x;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x + sqrt(lambda+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x - sqrt(lambda+n_x) * A.col(i);
  }

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.head(n_x) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x, n_x) = P;
  P_aug(5,5) = std_a * std_a;
  P_aug(6,6) = std_yawdd * std_yawdd;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  float scale_factor = sqrt(lambda + n_aug);
  Xsig_aug.col(0) = x_aug;
  
  for (int ii = 0; ii < n_aug; ii++){
      Xsig_aug.col(ii + 1) = x_aug + scale_factor * A.col(ii);
      Xsig_aug.col(ii + 1 + n_aug) = x_aug - scale_factor * A.col(ii);
  }

  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  
  
  //predict sigma points
  //avoid division by zero
  for (int ii = 0; ii < 2 * n_aug + 1; ii++){
      float px = Xsig_aug(0, ii);
      float py = Xsig_aug(1, ii);
      float v = Xsig_aug(2, ii);
      float psi = Xsig_aug(3, ii);
      float psi_dot = Xsig_aug(4, ii);
      float nu_a = Xsig_aug(5, ii);
      float nu_psi2 = Xsig_aug(6, ii);
      VectorXd x = VectorXd(5);
      x << px, py, v, psi, psi_dot;
      
      VectorXd term1 = VectorXd(5);
      VectorXd term2 = VectorXd(5);
      if(psi_dot != 0){
          float k_vs = v/psi_dot;
          float inside = psi + psi_dot * delta_t;
          term1(0) = k_vs * (sin(inside) - sin(psi));
          term1(1) = k_vs * (-cos(inside) + cos(psi));
          term1(2) = 0;
          term1(3) = psi_dot * delta_t;
          term1(4) = 0;
      }
      else{
          term1(0) = v * cos(psi) * delta_t;
          term1(1) = v * sin(psi) * delta_t;
          term1(2) = 0;
          term1(3) = 0;
          term1(4) = 0;
      }
      
      float dt2 = 0.5 * delta_t * delta_t;
      term2(0) = dt2 * cos(psi) * nu_a;
      term2(1) = dt2 * sin(psi) * nu_a;
      term2(2) = delta_t * nu_a;
      term2(3) = dt2 * nu_psi2;
      term2(4) = delta_t * nu_psi2;
      
      Xsig_pred.col(ii) = x + term1 + term2;
       
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
