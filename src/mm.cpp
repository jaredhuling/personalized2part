
#include "utils.h"
#include "thresholds.h"
#include "calc_U.h"


//[[Rcpp::export]]
Rcpp::List mmbcd_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
                     const Eigen::Map<Eigen::VectorXd> & Y,
                     const Eigen::Map<Eigen::VectorXi> & groups,
                     const Eigen::Map<Eigen::VectorXi> & unique_groups,
                     Eigen::VectorXd & group_weights,
                     Eigen::VectorXd & weights,
                     const Eigen::Map<Eigen::VectorXd> & offset,
                     Eigen::VectorXd & lambda,
                     const int &nlambda,
                     const double &lambda_min_ratio,
                     const double &alpha,
                     const double &tau,
                     const int &maxit,
                     const double &tol,
                     const bool &intercept,
                     std::vector<std::string> &family,
                     std::vector<std::string> &penalty)
{


    int ngroups = unique_groups.size();

    int nvars = X.cols();
    int nobs  = X.rows();

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


    double H_factor = 1.0;

    if (family[0] == "binomial")
    {
        H_factor = 0.25;
    }


    // set up threshold function
    thresh_func_ptr thresh_func;
    thresh_func = set_threshold_func(penalty[0]);


    // set up U function
    U_func_ptr U_func;
    U_func = set_U_func(family[0]);

    U_intercept_func_ptr U_intercept_func;
    U_intercept_func = set_U_intercept_func(family[0]);

    double gamma = 0.0;


    double b0 = 0.0;

    if (intercept)
    {
        if (family[0] == "gaussian")
        {
            b0 = (weights.array() * Y.array()).matrix().sum() / weights.sum();
        } else if (family[0] == "binomial")
        {
            double ybar = ( ((1.0 + Y.array()) * 0.5) * weights.array() ).matrix().sum() / weights.sum();
            b0 = log(ybar / (1.0 - ybar));
        } else if (family[0] == "gamma")
        {
            double ybar = (weights.array() * Y.array()).matrix().sum() / weights.sum();
            b0 = log(ybar);
        }
    }

    double b0_old = b0;




    // calculate X^TWX within each group
    std::vector<MatrixXd > x_list   = make_x_list(X, groups, ngroups, grp_idx);
    std::vector<MatrixXd > xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, weights);

    // compute largest eigenvalues within each group
    VectorXd eigenvals = compute_eigs(xtx_list);

    // set up default lambda sequence if no sequence provided
    if (nlambda_vec < 1)
    {
        lambda = setup_lambda(X, x_list, Y, weights, groups,
                              grp_idx, group_weights,
                              b0, U_func,
                              nlambda, lambda_min_ratio, penalty[0], alpha);
    }
    // END - set up default lambda

    nlambda_vec = lambda.size();

    VectorXi niter = VectorXi::Constant(nlambda_vec, maxit);




    VectorXd xbeta_cur(nobs);

    if (intercept)
    {
        xbeta_cur.array() = b0 + offset.array();
    } else
    {
        //xbeta_cur.setZero();
        xbeta_cur = offset;
    }

    VectorXd beta(nvars);
    VectorXd beta_old(nvars);

    MatrixXd beta_mat(nvars+1, nlambda_vec);

    beta.setZero();
    beta_old.setZero();
    beta_mat.setZero();


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

                b0_old = b0;

                double Upb = U_intercept_func(Y, weights, xbeta_cur, nobs) + b0;

                double total_penalty = 0.0;
                b0 = soft_thresh(Upb, total_penalty);

                xbeta_cur.array() += (b0 - b0_old);
            }
            beta_old = beta;

            if (family[0] == "gamma")
            {
                H_factor = (Y.array() / xbeta_cur.array().exp()).maxCoeff();
            }

            for (int g = 0; g < ngroups; ++g)
            {

                double stepsize = H_factor * eigenvals(g);

                int gr_size = grp_idx[g].size();
                VectorXd beta_subs(gr_size);
                for (int k = 0; k < gr_size; ++k)
                {
                    beta_subs(k) = beta(grp_idx[g][k]);
                }

                // calculate U vector
                VectorXd U_plus_beta = U_func(x_list[g], Y, weights, xbeta_cur, nobs).array() +
                    stepsize * beta_subs.array();


                //VectorXd beta_new = thresh_func(U_plus_beta, l1, gamma, l2, stepsize);


                double l1 = group_weights(g) * lam * tau;
                double lgr = group_weights(g) * lam * (1.0 - tau);

                VectorXd beta_new(gr_size);

                if (tau > 0.0)
                {
                    VectorXd beta_tmp(gr_size);
                    for (int k = 0; k < gr_size; ++k)
                    {
                        beta_tmp(k) = soft_thresh(U_plus_beta(k), l1);
                    }

                    beta_new = thresh_func(beta_tmp, lgr, lgr, stepsize);
                } else
                {
                    beta_new = thresh_func(U_plus_beta, lgr, lgr, stepsize);
                }

                bool anychanged = false;
                for (int k = 0; k < gr_size; ++k)
                {
                    if (beta_subs(k) != beta_new(k))
                    {
                        anychanged = true;
                        beta(grp_idx[g][k]) = beta_new(k);
                    }
                }


                // update residual if any estimate changed
                if (anychanged)
                {
                    xbeta_cur += (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());
                }

            }

            if (converged(beta, beta_old, tol))
            {
                niter(l) = i + 1;
                break;
            }

        } // end BCD iter loop

        beta_mat.col(l).tail(nvars) = beta;
        if (intercept)
        {
            beta_mat(0, l) = b0;
        }
    }

    return List::create(Named("beta")       = beta_mat,
                        Named("niter")      = niter,
                        Named("lambda")     = lambda,
                        Named("eigenvals")  = eigenvals,
                        Named("family")     = family[0],
                        Named("penalty")    = penalty[0]);
}
