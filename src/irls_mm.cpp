
#include "utils.h"
#include "thresholds.h"
#include "calc_U.h"


//[[Rcpp::export]]
Rcpp::List irls_mmbcd_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
                              const Eigen::Map<Eigen::VectorXd> & Y,
                              const Eigen::Map<Eigen::VectorXi> & groups,
                              const Eigen::Map<Eigen::VectorXi> & unique_groups,
                              Eigen::VectorXd & group_weights,
                              Eigen::VectorXd & weights,
                              Eigen::VectorXd & lambda,
                              const int &nlambda,
                              const double &lambda_min_ratio,
                              const double &alpha,
                              const double &tau,
                              const int &maxit,
                              const double &tol,
                              const int &maxit_irls,
                              const double &tol_irls,
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

    // set up threshold function
    thresh_func_ptr thresh_func;
    thresh_func = set_threshold_func(penalty[0]);


    // set up U function
    U_func_ptr U_func;
    U_func = set_U_func(family[0]);

    U_intercept_func_ptr U_intercept_func;
    U_intercept_func = set_U_intercept_func(family[0]);

    double gamma = 0.0;


    //double b0 = Y.sum() / double(nobs);

    VectorXd YY(nobs);
    if (family[0] == "binomial")
    {
        YY = (Y.array() + 1) / 2;
    }

    double b0 = 0;
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

    //b0 = Y.sum() / double(nobs);
    double b0_old = b0;

    // calculate X^TWX within each group
    std::vector<MatrixXd > x_list   = make_x_list(X, groups, ngroups, grp_idx);
    std::vector<MatrixXd > xtx_list(ngroups); // = compute_xtx_list(x_list, groups, ngroups, grp_idx, weights);

    // compute largest eigenvalues within each group
    VectorXd eigenvals(ngroups); // = compute_eigs(xtx_list);

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


    VectorXi niter = VectorXi::Constant(nlambda_vec, maxit_irls);



    VectorXd beta(nvars);
    VectorXd beta_old(nvars);
    VectorXd beta_irls_old(nvars);

    MatrixXd beta_mat(nvars+1, nlambda_vec);

    beta.setZero();
    beta_old.setZero();
    beta_mat.setZero();


    VectorXd p(nobs);
    VectorXd xbeta_cur(nobs);
    VectorXd resid_cur(nobs);
    VectorXd W(nobs);

    VectorXd deviance_vec(nlambda_vec);

    double deviance, deviance_old, weights_sum;


    if (intercept)
    {
        xbeta_cur.array() = b0;
    } else
    {
        xbeta_cur.setZero();
    }

    // loop over lamba values
    for (int l = 0; l < nlambda_vec; ++l)
    {
        double lam = lambda(l);
        deviance = 1e30;

        // irls iters
        for (int ii = 0; ii < maxit_irls; ++ii)
        {
            beta_irls_old = beta;
            deviance_old  = deviance;

            // -------------------------------------------------- //
            //                update quad approx                  //


            if (family[0] == "binomial")
            {
                // calculate mean function
                p = 1.0 / (1.0 + (-1.0 * xbeta_cur.array()).exp() );

                //std::cout << "min prob" << p.minCoeff() << std::endl;


                // construct weights and multiply by user-specified weights
                //W = weights.array() * p.array() * (1.0 - p.array());
                VectorXd var_vec = p.array() * (1.0 - p.array());
                W = weights.array() * var_vec.array();

                // make sure no weights are too small
                for (int k = 0; k < nobs; ++k)
                {
                    if (W(k) < 1e-5)
                    {
                        W(k) = 1e-5;
                    }
                }

                // update deviance
                deviance = 0.0;
                for (int ii = 0; ii < nobs; ++ii)
                {
                    if (Y(ii) == 1)
                    {
                        if (p(ii) > 1e-5)
                        {
                            deviance -= std::log(p(ii));
                        } else
                        {
                            // don't divide by zero
                            deviance -= std::log(1e-5);
                        }

                    } else
                    {
                        if (p(ii) <= 1.0 - 1e-5)
                        {
                            deviance -= std::log((1.0 - p(ii)));
                        } else
                        {
                            // don't divide by zero
                            deviance -= std::log(1.0 - 1e-5);
                        }

                    }
                }

                deviance *= 2.0;


                // here we update the residuals and multiply by user-specified weights, which
                // will be multiplied by X. ie X'resid_cur = X'Wz, where z is the working response from IRLS
                //resid_cur = weights.array() * (YY.array() - p.array()); // + xbeta_cur.array() * W.array().sqrt();
                resid_cur = (YY.array() - p.array()) / var_vec.array(); // + xbeta_cur.array() * W.array().sqrt();


            } else if (family[0] == "gamma")
            {
                p = Y.array() * (-1.0 * xbeta_cur.array()).exp();
                W = weights.array() * p.array();

                resid_cur = (p.array() - 1.0) / p.array();

                //std::cout << "W " << W.head(5).transpose() << std::endl;

                // make sure no weights are too small or too big
                for (int k = 0; k < nobs; ++k)
                {
                    if (W(k) < 1e-5)
                    {
                        W(k) = 1e-5;
                    } else if (W(k) > 1e5)
                    {
                        W(k) = 1e5;
                    }
                }

                deviance = (Y.array().log() - xbeta_cur.array() - p.array() + 1.0).matrix().sum();
                deviance *= -2.0;
            }

            if (converged_irls(deviance, deviance_old, tol_irls))
            {
                niter(l) = ii;
                break;
            }





            //Xsq = (W.array().sqrt().matrix().asDiagonal() * datX).array().square().colwise().sum();

            // this is needed for intercept updates
            weights_sum = W.sum();



            //weights_sum = weights.sum();

            //resid_cur = Y - xbeta_cur;


            //VectorXd realweights = W.array() * weights.array();
            xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, W);

            // compute largest eigenvalues within each group
            eigenvals = compute_eigs(xtx_list);

            // std::cout << "eigs " << eigenvals.head(5).transpose() << std::endl;

            //                                                    //
            // -------------------------------------------------- //

            /// bcd iters
            for (int i = 0; i < maxit; ++i)
            {
                // update intercept
                if (intercept)
                {

                    b0_old = b0;

                    double Upb = (W.array() * resid_cur.array() ).matrix().sum() / double(nobs) + b0; //  / ;

                    // double total_penalty = 0.0;
                    b0 = Upb; //soft_thresh(Upb, total_penalty);

                    //std::cout << "int U " << Upb << " int " << b0 << std::endl;

                    xbeta_cur.array() += (b0 - b0_old);
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

                    /*
                    VectorXd U_plus_beta = (x_list[g].transpose() *
                            ( W.array() * resid_cur.array()  ).matrix()).array() / double(nobs) +
                            eigenvals(g) * beta_subs.array();
                    */

                    VectorXd U_plus_beta(gr_size);
                    for (int k = 0; k < gr_size; ++k)
                    {
                        U_plus_beta(k) = ( X.col(grp_idx[g][k]).array() * W.array() * resid_cur.array()  ).matrix().sum() /
                            double(nobs) + eigenvals(g) * beta_subs(k);
                    }

                    double l1  = group_weights(g) * lam * tau;
                    double lgr = group_weights(g) * lam * (1.0 - tau);

                    VectorXd beta_new(gr_size);

                    if (tau > 0.0)
                    {
                        VectorXd beta_tmp(gr_size);
                        for (int k = 0; k < gr_size; ++k)
                        {
                            beta_tmp(k) = soft_thresh(U_plus_beta(k), l1);
                        }

                        beta_new = thresh_func(beta_tmp, lgr, gamma, lgr, eigenvals(g));
                    } else
                    {
                        beta_new = thresh_func(U_plus_beta, lgr, gamma, lgr, eigenvals(g));
                    }




                    /*
                    if (g == 0)
                    {
                        std::cout << "U update " << U_plus_beta.transpose() << std::endl;
                        std::cout << "beta update " << beta_new.transpose() << std::endl;
                    }
                     */


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
                        //xbeta_cur += (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());
                        //resid_cur -= (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());

                        for (int k = 0; k < gr_size; ++k)
                        {
                            xbeta_cur.array() += (X.col(grp_idx[g][k]).array() * (beta_new(k) - beta_subs(k)));
                            resid_cur.array() -= (X.col(grp_idx[g][k]).array() * (beta_new(k) - beta_subs(k)));
                        }
                    }

                }

                if (converged(beta, beta_old, tol))
                {
                    break;
                }

            } // end BCD iter loop

            // if (converged_irls(deviance, deviance_old, tol_irls))
            // {
            //     niter(l) = ii + 1;
            //     break;
            // }
            if (converged(beta, beta_irls_old, tol))
            {
                niter(l) = ii + 1;
                break;
            }
        } // end irls loop

        deviance_vec(l) = deviance;
        beta_mat.col(l).tail(nvars) = beta;
        if (intercept)
        {
            beta_mat(0, l) = b0;
        }
    }




    return List::create(Named("beta")      = beta_mat,
                        Named("niter")     = niter,
                        Named("lambda")    = lambda,
                        Named("deviance")  = deviance_vec,
                        Named("eigenvals") = eigenvals,
                        Named("family")    = family[0],
                        Named("penalty")   = penalty[0]);
}
