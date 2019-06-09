#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

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
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
    H_laser_ <<1,0,0,0,
            0,1,0,0;


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  //std::cout<<"bc this is measurement pack"<<measurement_pack.raw_measurements_ <<"with time stamp"<< measurement_pack.timestamp_<<std::endl;
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    float prev_t =  measurement_pack.timestamp_;//can optimize
    float dt = 500000.0/ 1000000.0;
    
    MatrixXd F(4,4);
      F<< 1,0,dt,0,
          0,1,0,dt,
          0,0,1,0,
          0,0,0,1;
 std::cout<<"F is "<<F<<std::endl;   
        MatrixXd Q(4,4);
    
float noise_a = 9;
    Q<<pow(dt,4)/4*pow(noise_a,2),0,pow(dt,3)/2*pow(noise_a,2),0,
      0,pow(dt,4)/4*pow(noise_a,2),0,pow(dt,3)/2*pow(noise_a,2),
      pow(dt,3)/2*pow(noise_a,2),0,pow(dt,2)*pow(noise_a,2),0,
      0,pow(dt,3)/2*pow(noise_a,2),0,pow(dt,2)*pow(noise_a,2);
    std::cout<<"Q is "<<Q<<std::endl;
     
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      VectorXd y = measurement_pack.raw_measurements_ - H_laser_*ekf_.x_;
      MatrixXd P_ = MatrixXd::Identity(4,4)*100;
     
      MatrixXd S = H_laser_*P_*H_laser_.transpose() + R_laser_;
      MatrixXd K = P_*H_laser_.transpose()*S.inverse();
      ekf_.x_ = ekf_.x_+ K*y;

      int size = (K*H_laser_).rows();
      MatrixXd I = MatrixXd::Identity(size,size);
      P_ = (I-K*H_laser_)*P_;

      ekf_.Init(ekf_.x_,P_,F,H_laser_,R_laser_,Q);
 
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates

  } else {
    // TODO: Laser updates

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
