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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
      //@todo: initialize state & cov mat
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      //@todo: initialize state & cov mat
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
}

void UKF::UpdateRadar(const MeasurementPackage &meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}