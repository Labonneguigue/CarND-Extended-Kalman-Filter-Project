#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

#define DBZ 0.0001    // Number below which it is considered 0

Tools::Tools(int size)
: mOverallSquaredSum(Eigen::VectorXd(size))
, counter(0)
{
    
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
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if ( ( estimations.size() > 0 ) && ( ground_truth.size() > 0 ) &&
         ( estimations.size() == ground_truth.size() ) )
    {
        assert(mOverallSquaredSum.size() == estimations[0].size());
        
        for ( ; counter < estimations.size() ; counter++){
            // sqrt ( sum ( ( gt(i) - est(i) )^2 ) / n )
            if ( ( estimations[counter].size() != ground_truth[counter].size() ) ||
                ( estimations[counter].size() != mOverallSquaredSum.size() ) ){
                cout << "\nProblem ! \n\n";
                return rmse;
            }
            else {
                Eigen::VectorXd difference = estimations[counter] - ground_truth[counter];
                difference = pow( difference.array(), 2);
                // I keep the overall sum of squared differences
                mOverallSquaredSum += difference;
            }
        }
        // I compute the result once the latest estimations has been added to the rolling sum
        rmse = mOverallSquaredSum / counter;
        rmse = rmse.array().sqrt();
        
        return rmse;
    }
    else
    {
        return rmse;
    }
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
