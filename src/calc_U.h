
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// pointer we will set to one of the U functions
typedef VectorXd (*U_func_ptr)(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> &y, VectorXd &weights, VectorXd &xbeta, int & n);

// pointer we will set to one of the U functions for intercept updating
typedef double (*U_intercept_func_ptr)(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);


U_func_ptr set_U_func(std::string & family);
U_intercept_func_ptr set_U_intercept_func(std::string & family);


// GAUSSIAN
VectorXd U_func_gaussian(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

double U_intercept_func_gaussian(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);


// BINOMIAL
VectorXd U_func_binomial(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

double U_intercept_func_binomial(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

// GAMMA
VectorXd U_func_gamma(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

double U_intercept_func_gamma(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

