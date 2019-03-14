
#include "fit_gamma.h"


//[[Rcpp::export]]
Rcpp::List fit_gamma_cpp(const Rcpp::NumericMatrix &X_,
                           const Rcpp::NumericVector &Y_,
                           const Rcpp::IntegerVector &groups_,
                           const Rcpp::IntegerVector &unique_groups_,
                           const Rcpp::NumericVector &group_weights_,
                           const Rcpp::NumericVector &weights_,
                           const Rcpp::NumericVector &offset_,
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
                           const bool strongrule)
{

    const MapMatd X  = Rcpp::as<MapMatd>(X_);
    const MapVecd Y  = Rcpp::as<MapVecd>(Y_);

    const MapVeci groups  = Rcpp::as<MapVeci>(groups_);
    const MapVeci unique_groups  = Rcpp::as<MapVeci>(unique_groups_);

    const MapVecd group_weights  = Rcpp::as<MapVecd>(group_weights_);
    const MapVecd weights  = Rcpp::as<MapVecd>(weights_);
    const MapVecd offset  = Rcpp::as<MapVecd>(offset_);
    const MapVecd lambda  = Rcpp::as<MapVecd>(lambda_);

    // set up object containing key parameters
    params pars;
    pars.tau = tau;
    pars.maxit = maxit;
    pars.maxit_irls = maxit_irls;
    pars.tol = tol;
    pars.tol_irls = tol_irls;
    pars.intercept_z = intercept;
    pars.intercept_s = intercept;
    pars.penalty = penalty[0];
    pars.nlambda = nlambda;
    pars.lambda_min_ratio = lambda_min_ratio;
    pars.strongrule = strongrule;


    gammamod gamma_obj(X, Y,
                       weights, offset,
                       lambda, groups, unique_groups,
                       group_weights,
                       pars);

    // initialization steop
    gamma_obj.initialize();

    // compute solution path
    VectorXi niter = gamma_obj.fit_path();


    // collect results
    MatrixXd beta             = gamma_obj.get_beta();
    VectorXd lamused          = gamma_obj.get_lambda();
    VectorXd deviance_vec     = gamma_obj.get_dev();
    VectorXd eigenvals        = gamma_obj.get_eigs();


    return List::create(Named("beta")      = beta,
                        Named("niter")     = niter,
                        Named("lambda")    = lamused,
                        Named("tau")       = tau,
                        Named("deviance")  = deviance_vec,
                        Named("eigenvals") = eigenvals,
                        Named("penalty")   = penalty[0]);
}
