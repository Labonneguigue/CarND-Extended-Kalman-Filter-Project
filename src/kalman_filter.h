#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:
    /**
     * Constructor
     */
    KalmanFilter();
    
    /**
     * Destructor
     */
    virtual ~KalmanFilter();
    
    /**
     * Init Initializes Kalman filter
     * @param x_in Initial state
     * @param P_in Initial state covariance
     * @param F_in Transition matrix
     * @param Q_in Process covariance matrix
     */
    void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in, Eigen::MatrixXd &Q_in);
    
    /**
     * Prediction Predicts the state and the state covariance
     * using the process model
     * @param dt Time between k and k+1 in s
     * @param noise_ax Model noise along x
     * @param noise_ay Model noise along y
     */
    void Predict(float dt, float noise_ax, float noise_ay);
    
    /**
     * Updates the state by using standard Kalman Filter equations
     * @param z The measurement at k+1
     */
    void Update(const Eigen::VectorXd &z, Eigen::MatrixXd H, Eigen::MatrixXd R);
    
    /**
     * Updates the state by using Extended Kalman Filter equations
     * @param z The measurement at k+1
     */
    void UpdateEKF(const Eigen::VectorXd &z, Eigen::MatrixXd H, Eigen::MatrixXd R);
    
    /**
     * Performs the common part of the update step.
     * @param y The error between the model and the new measurement.
     */
    void CommonUpdate(const Eigen::VectorXd &y);
    
    /**
     * Ensure that the angle is within the ]-PI;PI] range
     */
    void NormalizeAngle(double &angle);
    
    /**
     * Returns the state vector x
     */
    Eigen::VectorXd GetStateVector();
    
    /**
     * Returns the covariance matrix
     */
    Eigen::MatrixXd GetCovarianceMatrix();
    
private:
    
    bool isInitialized;
    
    // state vector
    Eigen::VectorXd x_;
    
    // state covariance matrix
    Eigen::MatrixXd P_;
    
    // state transition matrix
    Eigen::MatrixXd F_;
    
    // process covariance matrix
    Eigen::MatrixXd Q_;
    
    // measurement matrix
    Eigen::MatrixXd H_;
    
    // measurement covariance matrix
    Eigen::MatrixXd R_;
};

#endif /* KALMAN_FILTER_H_ */
