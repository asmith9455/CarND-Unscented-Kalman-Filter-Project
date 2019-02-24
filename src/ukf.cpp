#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = ::Eigen::VectorXd(2 * n_aug_ + 1);

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.rows(); ++i)
  {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  R_radar_ = Eigen::MatrixXd::Zero(3,3);
  R_laser_ = Eigen::MatrixXd::Zero(2,2);

  using ::std::pow;

  R_radar_(0,0) = pow(std_radr_, 2);
  R_radar_(1,1) = pow(std_radphi_, 2);
  R_radar_(2,2) = pow(std_radrd_, 2);

  R_laser_(0,0) = pow(std_laspx_, 2);
  R_laser_(1,1) = pow(std_laspy_, 2);

  /**
   * End DO NOT MODIFY section for measurement noise values 
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  std::cout << "\n\n-------------------------------" << std::endl;
  std::cout << "-------------------------------" << std::endl;
  std::cout << "PROCESS MEASUREMENT START" << std::endl;


  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_)
  {
    std::cout << "received a laser measurement, but settings say we shouldn't use it." << std::endl;
    return;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_)
  {
    std::cout << "received a radar measurement, but settings say we shouldn't use it." << std::endl;
    return;
  }

  if (!is_initialized_)
  {
    std::cout << "initializing" << std::endl;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      const double px = meas_package.raw_measurements_(0);
      const double py = meas_package.raw_measurements_(1);

      x_.fill(0.0);
      x_(0) = px;
      x_(1) = py;
      x_(3) = ::std::atan2(py, px);

      P_(0,0) = R_laser_(0,0);
      P_(1,1) = R_laser_(1,1);
      P_(2,2) = 1.0;
      P_(3,3) = 1.0;
      P_(4,4) = 1.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      const double rho = meas_package.raw_measurements_(0);
      const double phi = meas_package.raw_measurements_(1);
      //const double phidot = meas_package.raw_measurements_(2); unusable for initialization?

      x_(0) = rho * ::std::cos(phi);
      x_(1) = rho * ::std::sin(phi);
      x_(3) = 0.0;
      x_(3) = 0.0;
      x_(4) = 0.0;

      P_(0,0) = R_laser_(0,0);
      P_(1,1) = R_laser_(1,1);
      P_(2,2) = 1.0;
      P_(3,3) = 1.0;
      P_(4,4) = 1.0;
    }
    else
    {
      throw std::runtime_error("Unknown measurement type supplied");
    }
    is_initialized_ = true;

    time_us_ = meas_package.timestamp_;

    return;
  }

  const long long us_diff = meas_package.timestamp_ - time_us_;

  const double time_diff = static_cast<double>(us_diff) * 1e-6;

  std::cout << "time diff: " << time_diff << std::endl;

  Prediction(time_diff);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    throw std::runtime_error("Unknown measurement type supplied");
  }

  std::cout << "FINISHED PROCESSING MEASUREMENT" << std::endl;

  time_us_ = meas_package.timestamp_;

  std::cout << "new state: " << x_ << std::endl;
  std::cout << "new state cov: " << P_ << std::endl;
}

void UKF::Prediction(double delta_t)
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //generate sigma points

  // create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  /**
   * Student part begin
   */

  std::cout << "-------------------------------" << std::endl;
  std::cout << "PREDICTION START" << std::endl;

  std::cout << "generating augmented state" << std::endl;

  // create augmented mean state
  x_aug.block(0, 0, n_x_, 1) = x_;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.block(0, 0, n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd P_sqrt = P_aug.llt().matrixL();
  P_sqrt *= std::sqrt(lambda_ + n_aug_);

  // create augmented sigma points

  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, 2 * n_aug_).colwise() = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) += P_sqrt;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) -= P_sqrt;

  std::cout << "finished generating sigma points" << std::endl;
  std::cout << Xsig_aug << std::endl;

  //perform prediction

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  std::cout << "finished predicting sigma points" << std::endl;


  std::cout << "predicted sigma points: " << std::endl;

  std::cout << Xsig_pred_ << std::endl;

  x_.fill(0.0);
  P_.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    P_ += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }

  std::cout << "prediction complete" << std::endl;
}

void UKF::UpdateLidar(const MeasurementPackage &meas_package)
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  MatrixXd Zsig = MatrixXd::Zero(2, 2 * n_aug_ + 1);

  for (int i = 0; i < Xsig_pred_.cols(); ++i)
  {
      const double px = Xsig_pred_(0,i);
      const double py = Xsig_pred_(1,i);
      // const double v = Xsig_pred_(2,i);
      // const double psi = Xsig_pred_(3,i);
      // const double psi_dot = Xsig_pred_(4,i); //unused for radar update
      
      Zsig(0,i) = px;
      Zsig(1,i) = py; //@todo: verify assumption
      //@todo: make sure px, py are not less than eps (equal to 0)
  }

  VectorXd z_pred = VectorXd::Zero(2);
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      z_pred += Zsig.col(i) * weights_(i);
  }
  
  // calculate innovation covariance matrix S
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(2,2);
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      const Eigen::MatrixXd diff = Zsig.col(i) - z_pred;
      S += weights_(i) * diff * diff.transpose();
  }
  
  S += R_laser_;

  MatrixXd Tc = MatrixXd::Zero(n_x_, 2);

  if (Xsig_pred_.cols() != Zsig.cols())
  {
    std::runtime_error("sigma point count mismatch");
  }
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  // calculate Kalman gain K;
  
  Eigen::MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  
  x_ += K * (meas_package.raw_measurements_ - z_pred);
  
  P_ -= K * S * K.transpose();

}

void UKF::UpdateRadar(const MeasurementPackage &meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  MatrixXd Zsig = MatrixXd::Zero(3, 2 * n_aug_ + 1);

  for (int i = 0; i < Xsig_pred_.cols(); ++i)
  {
      const double px = Xsig_pred_(0,i);
      const double py = Xsig_pred_(1,i);
      const double v = Xsig_pred_(2,i);
      const double psi = Xsig_pred_(3,i);
      // const double psi_dot = Xsig_pred_(4,i); //unused for radar update
      
      using ::std::sqrt;
      using ::std::pow;
      using ::std::atan2;
      using ::std::cos;
      using ::std::sin;
      
      Zsig(0,i) = sqrt(pow(px,2) + pow(py,2));
      Zsig(1,i) = atan2(py, px); //@todo: verify assumption
      Zsig(2,i) = (px * cos(psi) * v + py * sin(psi) * v) / sqrt(pow(px,2) + pow(py,2));
      //@todo: make sure px, py are not less than eps (equal to 0)
  }

  VectorXd z_pred = VectorXd::Zero(3);
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      z_pred += Zsig.col(i) * weights_(i);
  }
  
  // calculate innovation covariance matrix S
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(3,3);
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      const Eigen::MatrixXd diff = Zsig.col(i) - z_pred;
      S += weights_(i) * diff * diff.transpose();
  }
  
  S += R_radar_;

  MatrixXd Tc = MatrixXd::Zero(n_x_, 3);

  if (Xsig_pred_.cols() != Zsig.cols())
  {
    std::runtime_error("sigma point count mismatch");
  }
  
  for(int i = 0; i < Zsig.cols(); ++i)
  {
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  // calculate Kalman gain K;
  
  Eigen::MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  
  x_ += K * (meas_package.raw_measurements_ - z_pred);
  
  P_ -= K * S * K.transpose();

}