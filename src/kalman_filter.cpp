#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  //std::cout<<"prdcit was called"<<std::endl;
  //std::cout<<"This is F"<<F_<<std::endl;
   //std::cout<<"This is x_"<<x_<<std::endl;
  x_ = F_*x_;
 
  P_=F_*P_*F_.transpose() +Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  
  std::cout<<"Kalman update was called"<<std::endl;
  VectorXd y = z - H_*x_;
  while (y(1) > 3.1415926) {
    y(1) -= 3.1415926;
}
while (y(1) < -3.1415926) {
    y(1) += 3.1415926;
}
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;

  int size = (K*H_).rows();
  MatrixXd I = MatrixXd::Identity(size,size);
  P_ = (I-K*H_)*P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  Tools tools;

  if (z[0]>0.05 && z[1]>0.05){
  MatrixXd Hj = tools.CalculateJacobian(x_);
  std::cout<<"Hj is"<<Hj<<std::endl;
   std::cout<<"x_ is is"<<x_<<std::endl;
 
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
   
  VectorXd hx(3);
    
  hx[0] = sqrt(pow(px,2)+pow(py,2));
  hx[1] = atan2(double(py),double(px));
  hx[2] = (px*vx+py*vy)/sqrt(pow(px,2)+pow(py,2));
  std::cout<<"this is measured hx[1]"<<std::endl;
  std::cout<<hx[1]<<std::endl;
  VectorXd y = z - hx;
    
while (y(1) > 3.1415926) {
    y(1) -= 3.1415926;
}
while (y(1) < -3.1415926) {
    y(1) += 3.1415926;
}
  MatrixXd S = Hj*P_*Hj.transpose() + R_;
  MatrixXd K = P_*Hj.transpose()*S.inverse();
  x_ = x_ + K*y;

  int size = (K*Hj).rows();
  MatrixXd I = MatrixXd::Identity(size,size);
  P_ = (I-K*Hj)*P_;
  }
  
 
}
