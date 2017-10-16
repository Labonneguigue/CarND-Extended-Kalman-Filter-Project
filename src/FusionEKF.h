#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
    /**
     * Constructor.
     */
    FusionEKF();
    
    /**
     * Destructor.
     */
    virtual ~FusionEKF();
    
    /**
     * Run the whole flow of the Kalman Filter from here.
     */
    void ProcessMeasurement(const MeasurementPackage &measurement_pack);
    /**
     * Return the instance of the Ext. Kalman Filter
     */
    KalmanFilter GetExtendedKalmanFilter();
    
private:
    // check whether the tracking toolbox was initialized or not (first measurement)
    bool is_initialized_;
    
    // previous timestamp
    long long previous_timestamp_;
    
    Eigen::MatrixXd R_laser_;
    Eigen::MatrixXd R_radar_;
    Eigen::MatrixXd H_laser_;
    Eigen::MatrixXd Hj_;
    
    float noise_ax_ = 9;
    float noise_ay_ = 9;
    
    KalmanFilter ekf_; ///< Kalman Filter update and prediction math lives in here.
};

#endif /* FusionEKF_H_ */
