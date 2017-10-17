#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//#define _MY_RMSE_ 0
#define DBZ 0.0001    // Number below which it is considered 0

Tools::Tools(int size)
: mOverallSquaredSum(Eigen::VectorXd(size))
, counter(0)
{
    mOverallSquaredSum << 0, 0, 0, 0;
    
    //@TODO: Move to Unit Test !
    
    vector<Eigen::VectorXd> estimations;
    VectorXd estimate(4);
    
    estimate(0) = 0.0F;
    estimate(1) = 1.0F;
    estimate(2) = 0.0F;
    estimate(3) = 0.0F;
    
    estimations.push_back(estimate);
    
    vector<Eigen::VectorXd> gts;
    VectorXd gt(4);
    
    gt(0) = -1.0F;
    gt(1) = -1.0F;
    gt(2) = 1.0F;
    gt(3) = 0.0F;
    
    gts.push_back(gt);
    
    Eigen::VectorXd result(4);
    
    result = CalculateRMSE(estimations, gts);
    
    cout << "RMSE : " << result(0) << " " << result(1) << " " << result(2) << " " << result(3) << "\n";
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
#ifdef _MY_RMSE_
    VectorXd rmse(4);
    
    int vectorSize = estimations.size();
    
    if ( ( vectorSize == 0 ) ||
         ( vectorSize != ground_truth.size() ) ||
         ( estimations[vectorSize - 1].size() != ground_truth[vectorSize - 1].size() ) )
    {
        std::cout << " ERROR \n";
        rmse << 0,0,0,0;
        return rmse;
    }
    else
    {
        Eigen::VectorXd difference(4);
        difference = estimations[vectorSize -1] - ground_truth[vectorSize -1];
        
        difference = difference.array() * difference.array();
        // I keep the overall sum of squared differences
        mOverallSquaredSum += difference;
        
        // I compute the result once the latest estimations has been added to the rolling sum
        rmse = mOverallSquaredSum / vectorSize;
        rmse = rmse.array().sqrt();
        
        return rmse;
    }
#else
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    
    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    
    //calculate the mean
    rmse = rmse/estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    //return the result
    return rmse;
#endif
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    
    Eigen::MatrixXd Hj(3,4);
    
    // retrieve state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    // Compute intermediate values used several time afterwards
    float v1 = ( px * px ) + ( py * py );
    
    if ( abs(v1) < DBZ ){
        v1 = DBZ;
    }
    
    float v2 = sqrt( v1 );
    float v3 = ( v1 * v2 );
    
    
    if ( abs(v3) < DBZ ){
        v3 = DBZ;
    }
    
    float H_2_0 = ( py * ( vx * py - vy * px ) ) / v3;
    float H_2_1 = ( px * ( px * vy - py * vx ) ) / v3;
    
    // Computation of Hj, the Jacobian matrix
    Hj <<  ( px / v2 ), ( py / v2 ), 0, 0,
    (-py / v1 ), ( px / v1 ), 0, 0,
    H_2_0 , H_2_1, ( px / v2 ), ( py / v2 );
    
    return Hj;
}
