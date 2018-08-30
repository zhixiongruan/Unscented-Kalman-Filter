Skip to content
 
Search or jump to…

Pull requests
Issues
Marketplace
Explore
 @zhixiongruan Sign out
1
0 0 stevemg9/udacity_p7_unscented_kalman_filter
 Code  Issues 0  Pull requests 0  Projects 0  Wiki  Insights
udacity_p7_unscented_kalman_filter/src/ukf.cpp
ea87eb7  13 days ago
 Steve Giardinelli Push to GitHub
     
Executable File  474 lines (350 sloc)  13.4 KB
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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

  // Initialization Flag
  is_initialized_ = false;

  // Time in us
  time_us_ = 0.0;

  // Sigma point spreading param
  lambda_ = 3 - n_x_;

  // State variable dimension
  n_x_ = 5;

  // Augmented state variable dimension
  n_aug_ = 7;

  // Initialize Matrix for Sigma Points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initialize Vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // If we have a RADAR or LiDAR measurement and the corresponding "use" flags are true:
  if((use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
    (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)){

    //--------------------------------------------------------------------------
    // If the UKF is not initialized, call the initialization function and skip 
    // update and predict steps
    //--------------------------------------------------------------------------
    if(!is_initialized_){
      InitUKF(meas_package);
      return;
    }

    //----------------------------------------------------------------------
    // If the UKF if already initialized, call the Predict and Update steps
    //----------------------------------------------------------------------

    // Calculate time elapsed since last measurement in seconds
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

    // Update the measurement time
    time_us_ = meas_package.timestamp_;

    // Call UKF Prediction Step
    Prediction(delta_t);

    // Call UKF Update Step for appropriate sensor
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Initializes the Unscented Kalman Filter by setting some of the state variables.
 * Function should be called on first measurement
 */
void UKF::InitUKF(MeasurementPackage meas_package){
  //Initialize state vector (x_) with default values for first measurement
  x_ << 1,
        1,
        1,
        1,
        0.1;

  // Initialize covariance matrix (P_) with default values
  P_ <<  0.15,     0,  0,  0,  0,
            0,  0.15,  0,  0,  0,
            0,     0,  1,  0,  0,
            0,     0,  0,  1,  0,
            0,     0,  0,  0,  1;

  // Initialize time with current timestamp
  time_us_ = meas_package.timestamp_;

  if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR){

  // Read in initial Radar measurements
  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];

  // Convert rho and phito px, and py.
  double px = rho * cos(phi);
  double py = rho * sin(phi);

  //Set Initial State
  x_(0) = px;
  x_(1) = py;
  }

  if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){

    x_(0) = meas_package.raw_measurements_[0];
    x_(1) = meas_package.raw_measurements_[1];
  }

  // Return from function after initialization
  is_initialized_ = true;
  return;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // ----------------------------------------------------------------------------
  // Sigma Point Generation
  //-----------------------------------------------------------------------------

  //Set Lambda for sigma points
  lambda_ = 3 - n_x_;

  // Declare Sigma Point Matrix
  MatrixXd X_sig = MatrixXd(n_x_, 2 * n_x_ + 1);

  // Calculate sqrt of Covariance Matrix (P_)
  MatrixXd P_sqrt = P_.llt().matrixL();

  // Calculate Sigma Points
  X_sig.col(0) = x_;

  for(int i = 0; i < n_x_; ++i){

    X_sig.col(i + 1)        = x_ + sqrt(lambda_ + n_x_) * P_sqrt.col(i);
    X_sig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * P_sqrt.col(i);
  }


  // ----------------------------------------------------------------------------
  // Sigma Point Augmentation
  //-----------------------------------------------------------------------------

  // ----- Declarations -----
  //Set lambda for n_aug_
  lambda_ = 3 - n_aug_;

  // Augmented Mean State Vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Augmented State Covariance Matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Augmented Sigma Point Matrix
  MatrixXd X_sig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Fill Augmented Mean State Vector
  x_aug.setZero();
  x_aug.head(5) = x_;

  // Fill Augmented State Covariance Matrix
  P_aug.setZero();
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //Calculate Square Root of P_aug
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  // Fill Augmented Sigma Point Matrix
  X_sig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; ++i){
    X_sig_aug.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * P_aug_sqrt.col(i);
    X_sig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * P_aug_sqrt.col(i);
  }


  // ----------------------------------------------------------------------------
  // Sigma Point Prediction
  //-----------------------------------------------------------------------------

  // 1/2 * dt^2
  double dt_sq = 0.5 * delta_t * delta_t;

  for(unsigned int i=0; i < X_sig_aug.cols(); ++i){

    // Get state variables from Augmented Sigma Point Matrix
    double px        = X_sig_aug(0,i);
    double py        = X_sig_aug(1,i);
    double v         = X_sig_aug(2,i);
    double yaw       = X_sig_aug(3,i);
    double yaw_d     = X_sig_aug(4,i);
    double nu_a      = X_sig_aug(5,i);
    double nu_yaw_dd = X_sig_aug(6,i);

    // Assign Vales to Predicted Sigma Point Matrix (Xsig_pred_)
    // Protect against division by zero
    if(abs(yaw_d) < 0.001){
      Xsig_pred_(0,i) = px + v * cos(yaw) * delta_t + dt_sq * cos(yaw) * nu_a;
      Xsig_pred_(1,i) = py + v * sin(yaw) * delta_t + dt_sq * sin(yaw) * nu_a;
    } else {
      Xsig_pred_(0,i) = px + (v / yaw_d) * (sin(yaw + yaw_d * delta_t) - sin(yaw))
                      + dt_sq * cos(yaw) * nu_a;

      Xsig_pred_(1,i) = py + (v / yaw_d) * (-cos(yaw + yaw_d * delta_t) + cos(yaw))
                      + dt_sq * sin(yaw) * nu_a;
    }

    Xsig_pred_(2,i) = v + delta_t * nu_a;
    Xsig_pred_(3,i) = yaw + yaw_d * delta_t + dt_sq * nu_yaw_dd;
    Xsig_pred_(4,i) = yaw_d + delta_t * nu_yaw_dd;
  }


  // ----------------------------------------------------------------------------
  // Predicted Mean and Covariance Assignment
  //-----------------------------------------------------------------------------

  // Set Weights
  for(int i = 0; i < 2 * n_aug_ + 1; ++i){
    if(i==0){
      weights_(i) = lambda_ / (lambda_ + double(n_aug_));
    } else {
      weights_(i) = 1.0 / (2.0 * (lambda_ + double(n_aug_)));
    }
  }

  // Predict Mean State Vector
  x_.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; ++i){

    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //Predict State Covariance Matrix
  P_.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; ++i){
    
    // State Vector Difference 
    VectorXd diff_x = Xsig_pred_.col(i) - x_;

    // Normalize Angle between -PI and PI
    while (diff_x(3)> M_PI) diff_x(3) -= 2.0*M_PI;
    while (diff_x(3)<-M_PI) diff_x(3) += 2.0*M_PI;

    P_ += weights_(i) * diff_x * diff_x.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Set dimensionality for LiDAR data
  int n_lidar = 2;

  // Extract LiDAR data into (z)
  VectorXd z = VectorXd(n_lidar);
  z = meas_package.raw_measurements_;

  // Declare Measurement Sigma Point Matrix (Z_sig)
  MatrixXd Z_sig = MatrixXd(n_lidar, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i){

    // Transform Predicted Sigma Points into the Measurement Space
    // Since Px and Py are measured directly no transform is required
    Z_sig(0,i) = Xsig_pred_(0,i);
    Z_sig(1,i) = Xsig_pred_(1,i);
  }

  // Calculate the Mean Predicted Measurement Vector (z_pred)
  VectorXd z_pred = VectorXd(n_lidar);
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i){
    z_pred += weights_(i) * Z_sig.col(i);
  }

  // Calculate the Innovation Covariance Matrix (S)
  MatrixXd R(n_lidar, n_lidar);
  R <<  std_laspx_ * std_laspx_,                        0,
                              0,  std_laspy_ * std_laspy_;

  MatrixXd S(n_lidar, n_lidar);
  S.setZero();

  for (int i = 0; i < 2 * n_aug_ + 1; ++i){

    VectorXd diff_z = Z_sig.col(i) - z_pred;

    S += weights_(i) * diff_z * diff_z.transpose();
  }

  S += R;

  // ------------UKF Update---------------

  // Declare Cross-correlation Matrix (Tc)
  MatrixXd Tc = MatrixXd(n_x_, n_lidar);

  // Calculate Cross-correlation Matrix
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i){
    
    VectorXd diff_X = Xsig_pred_.col(i) - x_;
    VectorXd diff_Z = Z_sig.col(i) - z_pred;

    Tc += weights_(i) * diff_X * diff_Z.transpose();
  }

  // Calculate the Kalman Gain (K)
  MatrixXd K = Tc * S.inverse();

  // Update State Mean Vector (x_) and Covariance Matrix (P_)
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // Set dimensionality for RADAR data
  int n_radar = 3;

  // Extract LiDAR data into (z)
  VectorXd z = VectorXd(n_radar);
  z = meas_package.raw_measurements_;

  // Declare Measurement Sigma Point Matrix (Z_sig)
  MatrixXd Z_sig = MatrixXd(n_radar, 2 * n_aug_ + 1);

  // Transform Predicted Sigma Points into the Measurement Space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i){

    //Extract State Variables from Predicted Sigma Point Matrix
    double px  = Xsig_pred_(0,i);
    double py  = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    // Calculate and assign values in measurement space
    double rho = sqrt(px*px + py*py);

    Z_sig(0,i) = rho;
    Z_sig(1,i) = atan2(py, px);
    Z_sig(2,i) = ((px * cos(yaw) * v) + (py * sin(yaw) * v)) / rho;
  }

  // Calculate the Mean Predicted Measurement Vector (z_pred)
  VectorXd z_pred = VectorXd(n_radar);
  z_pred.setZero();
  for (int i = 0; i <  2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Z_sig.col(i);
  }

  // Calculate the Innovation Covariance Matrix (S)
  MatrixXd R = MatrixXd(n_radar, n_radar);
  R <<  std_radr_ * std_radr_,                        0,                      0,
                            0,  std_radphi_*std_radphi_,                      0,
                            0,                        0,  std_radrd_*std_radrd_;

  MatrixXd S = MatrixXd(n_radar, n_radar);
  S.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd diff_z = Z_sig.col(i) - z_pred;

    while (diff_z(1)> M_PI) diff_z(1) -= 2.*M_PI;
    while (diff_z(1)<-M_PI) diff_z(1) += 2.*M_PI;

    S += weights_(i) * diff_z *diff_z.transpose();
  }

  S += R;

  // ------------UKF Update---------------

  // Declare Cross-correlation Matrix (Tc)
  MatrixXd Tc = MatrixXd(n_x_, n_radar);

  // Calculate Cross-correlation Matrix (Tc)
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i){
    
    VectorXd diff_X = Xsig_pred_.col(i) - x_;
    VectorXd diff_Z = Z_sig.col(i) - z_pred;

    // Normalize angles
    while (diff_X(3)> M_PI) diff_X(3) -= 2.*M_PI;
    while (diff_X(3)<-M_PI) diff_X(3) += 2.*M_PI;
    while (diff_Z(1)> M_PI) diff_Z(1) -= 2.*M_PI;
    while (diff_Z(1)<-M_PI) diff_Z(1) += 2.*M_PI;

    Tc += weights_(i) * diff_X * diff_Z.transpose();
  }

  // Calculate the Kalman Gain (K)
  MatrixXd K = Tc * S.inverse();

  // Calculate Difference to Normalize Angles
  VectorXd diff_z = z - z_pred;
  while (diff_z(1)> M_PI) diff_z(1) -= 2.*M_PI;
  while (diff_z(1)<-M_PI) diff_z(1) += 2.*M_PI;

  // Update State Mean Vector (x_) and Covariance Matrix (P_)
  x_ = x_ + K * diff_z;
  P_ = P_ - K * S * K.transpose();

}



© 2018 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
Press h to open a hovercard with more details.