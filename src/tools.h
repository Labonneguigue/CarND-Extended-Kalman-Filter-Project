#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
     * Constructor.
     */
    Tools(int size);
    
    /**
     * Destructor.
     */
    virtual ~Tools();
    
    /**
     * A helper method to calculate RMSE.
     */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                           const vector<VectorXd> &ground_truth);
    /**
     * Reset the counter and Vector keeping the data
     */
    void ResetRMSE();
    
    /**
     * A helper method to calculate Jacobians.
     */
    static MatrixXd CalculateJacobian(const VectorXd& x_state);
    
private:
    
    Eigen::VectorXd mOverallSquaredSum;
    int counter;
    
};

#endif /* TOOLS_H_ */
