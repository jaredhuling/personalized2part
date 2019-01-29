#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>

#include "calc_U.h"

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


std::vector<MatrixXd > compute_xtx_list(std::vector<MatrixXd > & x_list,
                                        const Eigen::Map<Eigen::VectorXi> & groups,
                                        const int &ngroups,
                                        std::vector<std::vector<int> > & grp_idx,
                                        Eigen::VectorXd & weights);

std::vector<MatrixXd > make_x_list(const Eigen::Map<Eigen::MatrixXd> & X,
                                   const Eigen::Map<Eigen::VectorXi> & groups,
                                   const int &ngroups,
                                   std::vector<std::vector<int> > & grp_idx);

Eigen::VectorXd compute_eigs(std::vector<MatrixXd > xtx_list);

Eigen::VectorXd setup_lambda(const Eigen::Map<Eigen::MatrixXd> & X,
                             std::vector<MatrixXd > & x_list,
                             const Eigen::Map<Eigen::VectorXd> & Y,
                             Eigen::VectorXd & weights,
                             const Eigen::Map<Eigen::VectorXi> & groups,
                             std::vector<std::vector<int> > & grp_idx,
                             Eigen::VectorXd & group_weights,
                             double b0,
                             U_func_ptr & U_func,
                             const int & nlambda,
                             const double & lambda_min_ratio,
                             std::string & penalty,
                             double & alpha);


Eigen::VectorXd setup_lambda(const Eigen::Map<Eigen::MatrixXd> & X,
                             const Eigen::Map<Eigen::VectorXd> & Y,
                             Eigen::VectorXd & weights,
                             const Eigen::Map<Eigen::VectorXi> & groups,
                             std::vector<std::vector<int> > & grp_idx,
                             Eigen::VectorXd & group_weights,
                             const int & nlambda,
                             const double & lambda_min_ratio,
                             std::string & penalty,
                             double & alpha);

std::vector<std::vector<int> > get_group_indexes(const VectorXi &groups,
                                                 const VectorXi &unique_groups,
                                                 const int &ngroups,
                                                 const int &nvars);

bool converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance);
