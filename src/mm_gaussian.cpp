
#include "utils.h"
#include "thresholds.h"
#include "calc_U.h"


//[[Rcpp::export]]
Rcpp::List mmbcd_gaussian_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
                              const Eigen::Map<Eigen::VectorXd> & Y,
                              const Eigen::Map<Eigen::VectorXi> & groups,
                              const Eigen::Map<Eigen::VectorXi> & unique_groups,
                              Eigen::VectorXd & group_weights,
                              Eigen::VectorXd & weights,
                              Eigen::VectorXd & lambda,
                              const int &nlambda,
                              const double &lambda_min_ratio,
                              const int &maxit,
                              const double &tol,
                              const bool &intercept)
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

    // setup threshold function

    thresh_func_ptr thresh_func;

    std::string penalty_val = "grp.lasso";

    thresh_func = set_threshold_func(penalty_val);

    double gamma = 0.0;
    double alpha = 1.0;

    // set up default lambda sequence if no sequence provided
    if (nlambda_vec < 1)
    {
        lambda = setup_lambda(X, Y, weights, groups, grp_idx, group_weights, nlambda, lambda_min_ratio,
                              penalty_val, alpha);
    }
    // END - set up default lambda

    nlambda_vec = lambda.size();

    // calculate X^TWX within each group
    std::vector<MatrixXd > x_list   = make_x_list(X, groups, ngroups, grp_idx);
    std::vector<MatrixXd > xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, weights);

    // compute largest eigenvalues within each group
    VectorXd eigenvals = compute_eigs(xtx_list);



    VectorXi niter = VectorXi::Constant(nlambda_vec, maxit);


    double b0 = Y.sum() / double(nobs);
    double b0_old = b0;

    VectorXd resid_cur(nobs);

    if (intercept)
    {
        resid_cur = Y.array() - b0;
    } else
    {
        resid_cur = Y;
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

                double Upb = (weights.array() * resid_cur.array()).matrix().sum() / double(nobs) + b0;

                double total_penalty = 0.0;
                b0 = soft_thresh(Upb, total_penalty);

                resid_cur.array() -= (b0 - b0_old);
            }
            beta_old = beta;

            for (int g = 0; g < ngroups; ++g)
            {
                int gr_size = grp_idx[g].size();
                VectorXd beta_subs(gr_size);
                for (int k = 0; k < gr_size; ++k)
                {
                    beta_subs(k) = beta(grp_idx[g][k]);
                }

                VectorXd U_plus_beta = (x_list[g].transpose() *
                        (weights.array() * resid_cur.array()).matrix()).array() / double(nobs) +
                        eigenvals(g) * beta_subs.array();

                double l1 = group_weights(g) * lam * alpha;
                double l2 = group_weights(g) * lam * (1.0 - alpha);

                VectorXd beta_new = thresh_func(U_plus_beta, l1, gamma, l2, eigenvals(g));

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
                    resid_cur -= (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());
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




    return List::create(Named("beta")      = beta_mat,
                        Named("niter")     = niter,
                        Named("lambda")    = lambda,
                        Named("eigenvals") = eigenvals);
}
