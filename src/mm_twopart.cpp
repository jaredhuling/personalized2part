#include "utils.h"
#include "thresholds_twopart.h"
#include "calc_U_twopart.h"


//[[Rcpp::export]]
Rcpp::List mmbcd_twopart_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
                             const Eigen::Map<Eigen::VectorXd> & Z,
                             const Eigen::Map<Eigen::MatrixXd> & Xs,
                             const Eigen::Map<Eigen::VectorXd> & S,
                             const Eigen::Map<Eigen::VectorXi> & groups,
                             const Eigen::Map<Eigen::VectorXi> & unique_groups,
                             Eigen::VectorXd & group_weights,
                             Eigen::VectorXd & weights,
                             Eigen::VectorXd & weights_s,
                             Eigen::VectorXd & lambda,
                             const int &nlambda,
                             const double &lambda_min_ratio,
                             const int &maxit,
                             const double &tol,
                             const bool &intercept,
                             std::vector<std::string> &penalty)
{


    int ngroups = unique_groups.size();

    int nvars = X.cols();
    int nobs  = X.rows();

    int nobs_s  = Xs.rows();

    int nlambda_vec = lambda.size();

    bool default_group_weights = bool(group_weights.size() < 1);

    // set up groups
    std::vector<std::vector<int> > grp_idx = get_group_indexes(groups, unique_groups, ngroups, nvars);

    if (default_group_weights)
    {
        group_weights.resize(ngroups);
        for (int g = 0; g < ngroups; ++g)
        {
            group_weights(g) = std::sqrt(double(grp_idx[g].size()));
        }
    }

    // END - set up groups


    double H_factor_binom = 0.25;
    double H_factor = 1.0;


    // set up threshold function
    thresh_func_twopart_ptr thresh_func;
    thresh_func = set_twopart_threshold_func(penalty[0]);

    std::string family = "binomial";
    std::string family_s = "gamma";


    // set up U function
    U_tp_func_ptr U_func;
    U_func = set_U_tp_func();

    U_intercept_func_ptr U_intercept_func;
    U_intercept_func = set_U_intercept_func(family);

    double gamma = 0.0;
    double alpha = 1.0;


    // set up intercept
    double b0 = 0.0;
    double b0_s = 0.0;

    if (intercept)
    {
        double zbar = ( ((1.0 + Z.array()) * 0.5) * weights.array() ).matrix().sum() / weights.sum(); //double(nobs);
        b0 = log(zbar / (1.0 - zbar));

        double sbar = (S.array() * weights_s.array()).matrix().sum() / weights_s.sum(); //double(nobs_s);
        //double sbar = S.sum() / double(nobs_s);
        b0_s = log(sbar);
    }

    double b0_old = b0;
    double b0_s_old = b0_s;

    std::cout << "int Z " << b0 << " int S " << b0_s << std::endl;


    // calculate X^TWX within each group
    //std::vector<MatrixXd > x_list   = make_x_list(X, groups, ngroups, grp_idx);
    //std::vector<MatrixXd > xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, weights);

    // compute largest eigenvalues within each group
    VectorXd eigenvals = compute_eigs_twopart(X, Xs);

    VectorXd penalty_adjustment(2);

    // set up default lambda sequence if no sequence provided
    if (nlambda_vec < 1)
    {
        VectorXd lambda_and_adjust = setup_lambda(X, Xs, Z, S, weights, weights_s, group_weights,
                                                  b0, b0_s, U_func,
                                                  nlambda, lambda_min_ratio, penalty[0], alpha);

        penalty_adjustment = lambda_and_adjust.head(2);

        std::cout << penalty_adjustment.transpose() << std::endl;

        lambda             = lambda_and_adjust.tail(nlambda);

    }

    penalty_adjustment.array() = 1.15;

    // END - set up default lambda

    nlambda_vec = lambda.size();

    VectorXi niter = VectorXi::Constant(nlambda_vec, maxit);




    VectorXd xbeta_cur(nobs);
    VectorXd xbeta_s_cur(nobs_s);

    if (intercept)
    {
        xbeta_cur.array() = b0;
        xbeta_s_cur.array() = b0_s;
    } else
    {
        xbeta_cur.setZero();
        xbeta_s_cur.setZero();
    }

    VectorXd beta(nvars);
    VectorXd beta_old(nvars);

    MatrixXd beta_mat(nvars+1, nlambda_vec);

    beta.setZero();
    beta_old.setZero();
    beta_mat.setZero();

    VectorXd beta_s(nvars);
    VectorXd beta_s_old(nvars);

    MatrixXd beta_s_mat(nvars+1, nlambda_vec);

    beta_s.setZero();
    beta_s_old.setZero();
    beta_s_mat.setZero();


    // loop over lamba values
    for (int l = 0; l < nlambda_vec; ++l)
    {
        double lam = lambda(l);

        /// bcd iters
        for (int i = 0; i < maxit; ++i)
        {
            // update intercept
            if (intercept)
            {

                b0_old   = b0;
                b0_s_old = b0_s;

                double Upb   = U_intercept_func_binomial(Z, weights,   xbeta_cur,  nobs) + b0;
                double Upb_s = U_intercept_func_gamma(S, weights_s, xbeta_s_cur, nobs_s) + b0_s;

                // double total_penalty = 0.0;
                b0   = Upb;   //soft_thresh(Upb,   total_penalty);
                b0_s = Upb_s; //soft_thresh(Upb_s, total_penalty);

                xbeta_cur.array()   += (b0 - b0_old);
                xbeta_s_cur.array() += (b0_s - b0_s_old);
            }
            beta_old   = beta;
            beta_s_old = beta_s;

            H_factor = (S.array() / xbeta_s_cur.array().exp()).maxCoeff();

            if (H_factor_binom > H_factor)
            {
                H_factor = H_factor_binom;
            }

            for (int g = 0; g < ngroups; ++g)
            {

                double stepsize = H_factor * eigenvals(g);

                VectorXd beta_subs(2);
                beta_subs(0) = beta(g);
                beta_subs(1) = beta_s(g);

                // calculate U vector
                VectorXd U_plus_beta = U_func(X.col(g), Xs.col(g), Z, S, weights, weights_s,
                                              xbeta_cur, xbeta_s_cur, nobs, nobs_s).array() +
                    stepsize * beta_subs.array();

                double l1 = group_weights(g) * lam * alpha;
                double l2 = group_weights(g) * lam * (1.0 - alpha);

                VectorXd beta_new = thresh_func(U_plus_beta, penalty_adjustment, l1, gamma, l2, stepsize);

                bool anychanged = false;
                for (int k = 0; k < 2; ++k)
                {
                    if (beta_subs(k) != beta_new(k))
                    {
                        anychanged = true;
                        beta(g)    = beta_new(0);
                        beta_s(g)  = beta_new(1);
                    }
                }


                // update residual if any estimate changed
                if (anychanged)
                {
                    xbeta_cur   += (X.col(g).array()  * (beta_new(0) - beta_subs(0))).matrix();
                    xbeta_s_cur += (Xs.col(g).array() * (beta_new(1) - beta_subs(1))).matrix();
                }

            }

            if (converged(beta, beta_old, tol) &&
                converged(beta_s, beta_s_old, tol))
            {
                niter(l) = i + 1;
                break;
            }

        } // end BCD iter loop

        beta_mat.col(l).tail(nvars)   = beta;
        beta_s_mat.col(l).tail(nvars) = beta_s;
        if (intercept)
        {
            beta_mat(0, l)   = b0;
            beta_s_mat(0, l) = b0_s;
        }
    }

    return List::create(Named("beta_z")     = beta_mat,
                        Named("beta_s")     = beta_s_mat,
                        Named("niter")      = niter,
                        Named("lambda")     = lambda,
                        Named("penalty_adjustment") = penalty_adjustment,
                        Named("eigenvals")  = eigenvals,
                        Named("penalty")    = penalty[0]);
}
