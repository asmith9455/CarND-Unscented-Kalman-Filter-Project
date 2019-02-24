#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.empty())
  {
    throw std::runtime_error("No estimations passed for calculation of rmse - the number of estimations must match the number of ground truth values.");
  }

  if (estimations.size() != ground_truth.size())
  {
    throw std::runtime_error("Incorrect number of estimations passed for calculation of rmse - the number of estimations must match the number of ground truth values.");
  }

  for (size_t i=0; i < estimations.size(); ++i) 
  {
    rmse += (estimations[i] - ground_truth[i]).array().pow(2).matrix();
  }
  
  rmse /= static_cast<double>(estimations.size());
  
  rmse = rmse.array().sqrt().matrix();
  
  return rmse;
}