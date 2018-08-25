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
 static long long previous_time_stamp;
 double delta_t;
 if (!is_initialized_) {
   
   cout << "Initializing" << endl;
   x_ << 5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;
   P_ << 0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << 5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      x_ << 5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;
    }

    //update the previous time stamp
    previous_time_stamp = meas_package.timestamp_;
    

    //done initializing, no need to predict or update
    is_initialized_ = true;
    return;
 }

 delta_t = (meas_package.timestamp_ - previous_time_stamp)/1000000.0;
 //cout << "Delta t: " << delta_t << endl;
 previous_time_stamp = meas_package.timestamp_;
 //cout << "Previous t: " << previous_time_stamp << endl;
 UKF::Prediction(delta_t);

 if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      cout << "In Lidar" << endl;
      //UKF::UpdateLidar(meas_package);
 }
 else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
   cout << "In Radar" << endl;
      UKF::UpdateRadar(meas_package);
 }

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

  cout << "IN Predict" << endl;
  
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_x;

  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  
  //set first colume of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x_ - sqrt(lambda+n_x) * A.col(i);
  }
  //cout << Xsig << "\n" << endl;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.head(n_x) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x, n_x) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  //cout << P_aug << "\n" << endl;

    //create augmented sigma points
  float scale_factor = sqrt(lambda + n_aug);
  Xsig_aug.col(0) = x_aug;

  MatrixXd A_aug = P_aug.llt().matrixL();
  
  for (int ii = 0; ii < n_aug; ii++){
      Xsig_aug.col(ii + 1) = x_aug + scale_factor * A_aug.col(ii);
      Xsig_aug.col(ii + 1 + n_aug) = x_aug - scale_factor * A_aug.col(ii);
  }

 Xsig_pred_ = MatrixXd(n_x, 2 * n_aug + 1);
  
  
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
      
      Xsig_pred_.col(ii) = x + term1 + term2;
      //cout << Xsig_pred_ << "\n" << endl;
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
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;
  cout << "Lambda: " << lambda << endl;
  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    cout << "v2: " << v2 << endl;
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  
  S = S + R;
  cout << "S: " << S << endl;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
