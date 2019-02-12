
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// pointer we will set to one of the thresholding functions
typedef VectorXd (*thresh_func_ptr)(VectorXd &value, double &penalty, double &l2, double &denom);

thresh_func_ptr set_threshold_func(std::string & penalty);

VectorXd block_soft_thresh(VectorXd & a, double & penalty, double &l2, double &denom);

VectorXd coop_block_soft_thresh(VectorXd & a, double & penalty, double &l2, double &denom);

double soft_thresh(double & a, double & penalty);

VectorXd phi_j_v(VectorXd & v, int & j);
