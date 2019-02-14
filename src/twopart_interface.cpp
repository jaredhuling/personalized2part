
#include "twopart.h"


//[[Rcpp::export]]
Rcpp::List fit_twopart_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
                           const Eigen::Map<Eigen::VectorXd> & Z,
                           const Eigen::Map<Eigen::MatrixXd> & Xs,
                           const Eigen::Map<Eigen::VectorXd> & S,
                           const Eigen::Map<Eigen::VectorXi> & groups,
                           const Eigen::Map<Eigen::VectorXi> & unique_groups,
                           Eigen::VectorXd & group_weights,
                           Eigen::VectorXd & weights,
                           Eigen::VectorXd & weights_s,
                           const Eigen::Map<Eigen::VectorXd> & offset,
                           const Eigen::Map<Eigen::VectorXd> & offset_s,
                           Eigen::VectorXd & lambda,
                           const int &nlambda,
                           const double &lambda_min_ratio,
                           const double &tau,
                           const int &maxit,
                           const double &tol,
                           const int &maxit_irls,
                           const double &tol_irls,
                           const bool &intercept,
                           std::vector<std::string> &penalty,
                           const bool &opposite_signs)
{

    // set up object containing key parameters
    params pars;
    pars.tau = tau;
    pars.maxit = maxit;
    pars.maxit_irls = maxit_irls;
    pars.tol = tol;
    pars.tol_irls = tol_irls;
    pars.intercept = intercept;
    pars.penalty = penalty[0];
    pars.opposite_signs = opposite_signs;
    pars.nlambda = nlambda;
    pars.lambda_min_ratio = lambda_min_ratio;


    twopart tp_obj(X, Xs, Z, S,
                   weights, weights_s, offset, offset_s,
                   lambda, groups, unique_groups,
                   group_weights,
                   pars);

    // initialization steop
    tp_obj.initialize();

    // compute solution path
    VectorXi niter = tp_obj.fit_path();


    // collect results
    MatrixXd beta_z             = tp_obj.get_beta_z();
    MatrixXd beta_s             = tp_obj.get_beta_s();
    VectorXd lamused            = tp_obj.get_lambda();
    VectorXd deviance_z_vec     = tp_obj.get_dev_z();
    VectorXd deviance_s_vec     = tp_obj.get_dev_s();
    VectorXd penalty_adjustment = tp_obj.get_pen_adj();
    VectorXd eigenvals          = tp_obj.get_eigs();


    return List::create(Named("beta_z")      = beta_z,
                        Named("beta_s")      = beta_s,
                        Named("niter")       = niter,
                        Named("lambda")      = lamused,
                        Named("tau")         = tau,
                        Named("deviance_z")  = deviance_z_vec,
                        Named("deviance_s")  = deviance_s_vec,
                        Named("penalty_adjustment") = penalty_adjustment,
                        Named("eigenvals") = eigenvals,
                        Named("penalty")   = penalty[0]);
}
