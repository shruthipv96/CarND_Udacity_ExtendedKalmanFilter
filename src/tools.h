#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
   * A helper function to predict the Radar measurements in polar co-ordinates
   */
  Eigen::VectorXd PredictRadarMeasurement(double px,double py,double vx,double vy);

  /**
    * A helper function to convert from Polar to Cartesian Cordinate
    */
  Eigen::VectorXd ConvertPolarToCartesian(float ro,float theta,float ro_dot);


};

#endif  // TOOLS_H_
