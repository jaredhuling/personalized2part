
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// pointer we will set to one of the thresholding functions
typedef VectorXd (*thresh_func_twopart_ptr)(VectorXd &value, VectorXd & penalty_factor, double &penalty, double &l1, double &denom);

thresh_func_twopart_ptr set_twopart_threshold_func(std::string & penalty);

VectorXd block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & penalty, double &l1, double &denom);

VectorXd coop_block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & penalty, double &l1, double &denom);

