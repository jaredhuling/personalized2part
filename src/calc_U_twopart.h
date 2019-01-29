
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// pointer we will set to one of the U functions
typedef VectorXd (*U_tp_func_ptr)(const VectorXd &x_col, const VectorXd &x_col_s,
                  const Eigen::Map<Eigen::VectorXd> &z,
                  const Eigen::Map<Eigen::VectorXd> &s,
                  VectorXd &weights,
                  VectorXd &weights_s,
                  VectorXd &xbeta,
                  VectorXd &xbeta_s,
                  int & n,
                  int & n_s);

// pointer we will set to one of the U functions for intercept updating
typedef double (*U_tp_intercept_func_ptr)(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);


U_tp_func_ptr set_U_tp_func();
U_tp_intercept_func_ptr set_U_tp_intercept_func(std::string & family);


//double U_intercept_func_gaussian(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);


// TWO PART
VectorXd U_func_twopart(const VectorXd &x_col, const VectorXd &x_col_s,
                        const Eigen::Map<Eigen::VectorXd> &z,
                        const Eigen::Map<Eigen::VectorXd> &s,
                        VectorXd &weights,
                        VectorXd &weights_s,
                        VectorXd &xbeta,
                        VectorXd &xbeta_s,
                        int & n,
                        int & n_s);

//double U_intercept_func_binomial(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

// GAMMA
VectorXd U_func_gamma(const VectorXd &x_col, const VectorXd &x_col_s,
                      const Eigen::Map<Eigen::VectorXd> &z,
                      const Eigen::Map<Eigen::VectorXd> &s,
                      VectorXd &weights,
                      VectorXd &weights_s,
                      VectorXd &xbeta,
                      VectorXd &xbeta_s,
                      int & n,
                      int & n_s);

//double U_intercept_func_gamma(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n);

