#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    if ( ( estimations.size() > 0 ) && ( ground_truth.size() > 0 ) )
    {
    
        VectorXd result( estimations[0].size() );
        
        for ( int i=0 ; i < std::max(estimations.size(),ground_truth.size()) ; i++){
            
        }
        
        if (estimations.size() != ground_truth.size()){
            // If size is different, I return an empty Vector of the size
            // of the ground truth vector.
            
            //result << 0, 0, 0, 0;
        }
        
        return result;
    }
    else
    {
        Eigen::VectorXd result(1);
        result << -1;
        return result ;
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
    float v2 = sqrt( v1 );
    float v3 = ( v1 * v2 );
    
    // Check if division by zero is occuring
    if ( abs( v1 ) < 0.0001 ){
        cout << "Division by 0 in CalculateJacobian() function.\n";
        return Hj;
    }
    
    float H_2_0 = ( py * ( vx * py - vy * px ) ) / v3;
    float H_2_1 = ( px * ( px * vy - py * vx ) ) / v3;
    
    // Computation of Hj, the Jacobian matrix
    Hj <<  ( px / v2 ), ( py / v2 ), 0, 0,
           (-py / v1 ), ( px / v1 ), 0, 0,
           H_2_0 , H_2_1, ( px / v2 ), ( py / v2 );
    
    return Hj;
}
