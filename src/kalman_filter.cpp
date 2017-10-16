#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define DBZ 0.0001    // Number below which it is considered 0

KalmanFilter::KalmanFilter()
    : isInitialized(false)
    {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    Q_ = Q_in;
    isInitialized = true;
}

void KalmanFilter::Predict(float dt, float noise_ax, float noise_ay) {
    /**
     TODO:
     * predict the state
     */
    F_ <<  1, 0, dt, 0,
    0, 1, 0, dt,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
    Q_ << pow(dt, 4) * noise_ax / 4 , 0, pow(dt, 3) * noise_ax / 2, 0,
          0, pow(dt, 4) * noise_ay / 4 , 0, pow(dt, 3) * noise_ay / 2,
          pow(dt, 3) * noise_ax / 2 , 0, pow(dt, 2) * noise_ax, 0,
          0, pow(dt, 3) * noise_ay / 2, 0, pow(dt, 2) * noise_ay;
    
    // u is zero in our case.
    assert(isInitialized);
    x_ = F_ * x_ ;// + u;
    MatrixXd Ft = F_.transpose();
    P_ = ( F_ * P_ * Ft ) + Q_;
}

void KalmanFilter::Update(const VectorXd &z, Eigen::MatrixXd H, Eigen::MatrixXd R) {
    /**
     TODO:
     * update the state by using Kalman Filter equations
     */
    assert(isInitialized);
    H_ = H;
    R_ = R;
    VectorXd y = z - ( H_ * x_ );
    CommonUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, Eigen::MatrixXd H, Eigen::MatrixXd R) {
    /**
     TODO:
     * update the state by using Extended Kalman Filter equations
     */
    /*
     Since z is in polar coordinates, I need to transpose the
     x state vector from the cartesian coordiantes to the
     polar ones
     */
    assert(isInitialized);
    H_ = H;
    R_ = R;
    VectorXd h = Eigen::VectorXd(3); // h(x_)
    
    h(0) = sqrt( x_(0) * x_(0) + x_(1) * x_(1) ); // rho
    h(1) = atan2( x_(1), x_(0) ); // theta
    h(2) = ( x_(0) * x_(2) + x_(1) * x_(3) ) / h(0); // rho dot
    
    VectorXd y = z - h;
    NormalizeAngle( y(1) );
    CommonUpdate(y);
}

void KalmanFilter::CommonUpdate(const Eigen::VectorXd &y){
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    //new state
    x_ = x_ + (K * y);
    
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size); // Identity matrix
    
    P_ = (I - ( K * H_ ) ) * P_;
}

void KalmanFilter::NormalizeAngle(double &angle){
    while ( angle > M_PI ){
        angle -= ( 2 * M_PI );
    }
    while ( angle < -M_PI ){
        angle += ( 2 * M_PI );
    }
}

Eigen::VectorXd KalmanFilter::GetStateVector(){
    return x_;
}

Eigen::MatrixXd KalmanFilter::GetCovarianceMatrix(){
    return P_;
}
