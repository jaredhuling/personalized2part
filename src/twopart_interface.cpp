
#include "twopart.h"


//[[Rcpp::export]]
Rcpp::List fit_twopart_cpp(const Rcpp::NumericMatrix &X_,
                           const Rcpp::NumericVector &Z_,
                           const Rcpp::NumericMatrix &Xs_,
                           const Rcpp::NumericVector &S_,
                           const Rcpp::IntegerVector &groups_,
                           const Rcpp::IntegerVector &unique_groups_,
                           const Rcpp::NumericVector &group_weights_,
                           const Rcpp::NumericVector &weights_,
                           const Rcpp::NumericVector &weights_s_,
                           const Rcpp::NumericVector &offset_,
                           const Rcpp::NumericVector &offset_s_,
                           const Rcpp::NumericVector &lambda_,
                           const int nlambda,
                           const double lambda_min_ratio,
                           const double tau,
                           const int maxit,
                           const double tol,
                           const int maxit_irls,
                           const double tol_irls,
                           const bool intercept,
                           const std::vector<std::string> penalty,
                           const bool opposite_signs,
                           const bool strongrule)
{

    const MapMatd X  = Rcpp::as<MapMatd>(X_);
    const MapVecd Z  = Rcpp::as<MapVecd>(Z_);
    const MapMatd Xs = Rcpp::as<MapMatd>(Xs_);
    const MapVecd S  = Rcpp::as<MapVecd>(S_);

    const MapVeci groups  = Rcpp::as<MapVeci>(groups_);
    const MapVeci unique_groups  = Rcpp::as<MapVeci>(unique_groups_);

    const MapVecd group_weights  = Rcpp::as<MapVecd>(group_weights_);
    const MapVecd weights  = Rcpp::as<MapVecd>(weights_);
    const MapVecd weights_s  = Rcpp::as<MapVecd>(weights_s_);
    const MapVecd offset  = Rcpp::as<MapVecd>(offset_);
    const MapVecd offset_s  = Rcpp::as<MapVecd>(offset_s_);
    const MapVecd lambda  = Rcpp::as<MapVecd>(lambda_);

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
    pars.strongrule = strongrule;


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
    double scale_pos            = tp_obj.get_scale();


    return List::create(Named("beta_z")      = beta_z,
                        Named("beta_s")      = beta_s,
                        Named("niter")       = niter,
                        Named("lambda")      = lamused,
                        Named("tau")         = tau,
                        Named("deviance_z")  = deviance_z_vec,
                        Named("deviance_s")  = deviance_s_vec,
                        Named("penalty_adjustment") = penalty_adjustment,
                        Named("likelihood_scale_factor") = scale_pos,
                        Named("eigenvals") = eigenvals,
                        Named("penalty")   = penalty[0]);
}
