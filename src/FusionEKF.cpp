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
    H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
    
    Hj_ = MatrixXd(3, 4);
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0,      0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0,      0,
    0,    0.0009, 0,
    0,    0,      0.09;
    
    /**
     TODO:
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    
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
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        cout << "EKF: " << endl;
        VectorXd x = VectorXd(4);
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            assert( measurement_pack.raw_measurements_.size() == 3 );
            float rho = measurement_pack.raw_measurements_[0]; // Distance
            float phi = measurement_pack.raw_measurements_[1]; // Bearing
            x(0) = rho * cos(phi);
            x(1) = rho * sin(phi);
            x(2) = 0.0F;
            x(3) = 0.0F;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            assert( measurement_pack.raw_measurements_.size() == 2 );
            x(0) = measurement_pack.raw_measurements_[0];
            x(1) = measurement_pack.raw_measurements_[1];
            x(2) = 0.0F;
            x(3) = 0.0F;
        }
        
        // Initialization of the covariance matrix
        Eigen::MatrixXd P = Eigen::MatrixXd(4, 4);
        P << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
        
        // State transition matrix construction
        MatrixXd F = MatrixXd(4, 4);
        F << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
        
        // Complete process noise covariance matrix construction
        MatrixXd Q = MatrixXd(4, 4);
        Q << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
        
        // Initialization of the Extended Kalman Filter instance
        ekf_.Init(x, P, F, Q);
        
        
        /* Store the timestamp to be able to calculate dt when
         the next measurement comes.
         */
        previous_timestamp_ = measurement_pack.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

        /**
         TODO:
         * Update the state transition matrix F according to the new elapsed time.
         - Time is measured in seconds.
         * Update the process noise covariance matrix.
         * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
         */
        
        //compute the time elapsed between the current and previous measurements
        float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
        previous_timestamp_ = measurement_pack.timestamp_;
    
    if ( dt > 0.001){
    
        ekf_.Predict(dt, noise_ax_, noise_ay_);
        
    }
        /*****************************************************************************
         *  Update
         ****************************************************************************/
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /* Radar updates
             Since we have ( Distance, Bearing and Distance rate ) as x_ vector
             we need to map them to the ( px, py, vx, vy ) state vector.
             This introduces non-linearity. It is solved by linearization
             around the current position.
             */
            Eigen::MatrixXd H = Tools::CalculateJacobian(ekf_.GetStateVector());
            ekf_.UpdateEKF(measurement_pack.raw_measurements_, H, R_radar_);
            
        } else {
            // Laser updates
            ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
        }
        
        // print the output
        cout << "x_ = " << ekf_.GetStateVector() << endl;
        cout << "P_ = " << ekf_.GetCovarianceMatrix() << endl;
}

KalmanFilter FusionEKF::GetExtendedKalmanFilter(){
    return ekf_;
}
