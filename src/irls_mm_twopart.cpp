
#include "utils.h"
#include "thresholds_twopart.h"
#include "thresholds.h"
#include "calc_U_twopart.h"


//[[Rcpp::export]]
Rcpp::List irls_mmbcd_twopart_cpp(const Eigen::Map<Eigen::MatrixXd> & X,
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
                                  const double &tau,
                                  const int &maxit,
                                  const double &tol,
                                  const int &maxit_irls,
                                  const double &tol_irls,
                                  const bool &intercept,
                                  std::vector<std::string> &penalty,
                                  const bool &opposite_signs)
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


    //double b0 = Y.sum() / double(nobs);

    VectorXd ZZ = (Z.array() + 1) / 2;

    double b0   = 0.0;
    double b0_s = 0.0;
    if (intercept)
    {
        double zbar = ( (ZZ.array() * weights.array()) ).matrix().sum() /  weights.sum(); // double(nobs);
        b0 = log(zbar / (1.0 - zbar));

        double sbar = (S.array() * weights_s.array()).matrix().sum() / weights_s.sum(); // double(nobs_s);
        //double sbar = S.sum() / double(nobs_s);
        b0_s = log(sbar);
    }

    //b0 = Y.sum() / double(nobs);
    double b0_old = b0;
    double b0_s_old = b0_s;

    // calculate X^TWX within each group
    //std::vector<MatrixXd > x_list   = make_x_list(X, groups, ngroups, grp_idx);
    //std::vector<MatrixXd > xtx_list(ngroups); // = compute_xtx_list(x_list, groups, ngroups, grp_idx, weights);

    // compute largest eigenvalues within each group
    VectorXd eigenvals(ngroups); // = compute_eigs(xtx_list);



    // set up default lambda sequence if no sequence provided
    VectorXd penalty_adjustment(2);

    // set up default lambda sequence if no sequence provided
    if (nlambda_vec < 1)
    {
        VectorXd lambda_and_adjust = setup_lambda(X, Xs, Z, S, weights, weights_s, group_weights,
                                                  b0, b0_s, U_func,
                                                  nlambda, lambda_min_ratio, penalty[0], alpha);


        penalty_adjustment = lambda_and_adjust.head(2);

        // std::cout << penalty_adjustment.transpose() << std::endl;

        lambda             = lambda_and_adjust.tail(nlambda);

    }

    penalty_adjustment.array() = 1.01;
    // END - set up default lambda

    nlambda_vec = lambda.size();


    VectorXi niter = VectorXi::Constant(nlambda_vec, maxit_irls);



    VectorXd beta(nvars);
    VectorXd beta_old(nvars);
    VectorXd beta_irls_old(nvars);

    VectorXd beta_s(nvars);
    VectorXd beta_s_old(nvars);
    VectorXd beta_irls_s_old(nvars);

    MatrixXd beta_mat(nvars+1, nlambda_vec);
    MatrixXd beta_s_mat(nvars+1, nlambda_vec);

    beta.setZero();
    beta_old.setZero();
    beta_mat.setZero();

    beta_s.setZero();
    beta_s_old.setZero();
    beta_s_mat.setZero();


    VectorXd p(nobs);
    VectorXd p_s(nobs_s);
    VectorXd xbeta_cur(nobs);
    VectorXd resid_cur(nobs);
    VectorXd xbeta_s_cur(nobs_s);
    VectorXd resid_s_cur(nobs_s);
    VectorXd W(nobs);
    VectorXd W_s(nobs_s);

    VectorXd deviance_vec(nlambda_vec);
    VectorXd deviance_s_vec(nlambda_vec);

    double deviance, deviance_old, weights_sum, weights_s_sum, deviance_s, deviance_s_old;

    double mult_1 = 1.0;
    if (opposite_signs)
    {
        mult_1 = -1.0;
    }


    if (intercept)
    {
        xbeta_cur.array()   = b0;
        xbeta_s_cur.array() = b0_s;
    } else
    {
        xbeta_cur.setZero();
        xbeta_s_cur.setZero();
    }

    // loop over lamba values
    for (int l = 0; l < nlambda_vec; ++l)
    {
        double lam = lambda(l);
        deviance   = 1e30;
        deviance_s = 1e30;

        // irls iters
        for (int ii = 0; ii < maxit_irls; ++ii)
        {
            beta_irls_old   = beta;
            deviance_old    = deviance;

            beta_irls_s_old = beta_s;
            deviance_s_old  = deviance_s;

            // -------------------------------------------------- //
            //                update quad approx                  //


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
                if (Z(ii) == 1)
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
            resid_cur = (ZZ.array() - p.array()) / var_vec.array(); // + xbeta_cur.array() * W.array().sqrt();

            /////// ----------------
            ///////   GAMMA UPDATE
            /////// ----------------

            p_s = S.array() * (-1.0 * xbeta_s_cur.array()).exp();
            W_s = weights_s.array() * p_s.array();

            resid_s_cur = (p_s.array() - 1.0) / p_s.array();

            //std::cout << "W " << W.head(5).transpose() << std::endl;

            // make sure no weights are too small or too big
            for (int k = 0; k < nobs_s; ++k)
            {
                if (W_s(k) < 1e-5)
                {
                    W_s(k) = 1e-5;
                } else if (W_s(k) > 1e5)
                {
                    W_s(k) = 1e5;
                }
            }

            deviance_s = (S.array().log() - xbeta_s_cur.array() - p_s.array() + 1.0).matrix().sum();
            deviance_s *= -2.0;

            if (converged_irls(deviance,   deviance_old,   tol_irls) &&
                converged_irls(deviance_s, deviance_s_old, tol_irls))
            {
                niter(l) = ii;
                break;
            }





            //Xsq = (W.array().sqrt().matrix().asDiagonal() * datX).array().square().colwise().sum();

            // this is needed for intercept updates
            weights_sum   = W.sum();
            weights_s_sum = W_s.sum();


            //weights_sum = weights.sum();

            //resid_cur = Y - xbeta_cur;


            //VectorXd realweights = W.array() * weights.array();
            //xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, W);

            // compute largest eigenvalues within each group
            //eigenvals = compute_eigs(xtx_list);
            VectorXd eigenvals = compute_eigs_twopart(X, Xs, W, W_s);

            //std::cout << "eigs " << eigenvals.head(5).transpose() << std::endl;

            //                                                    //
            // -------------------------------------------------- //

            /// bcd iters
            for (int i = 0; i < maxit; ++i)
            {
                // update intercept
                if (intercept)
                {

                    b0_old   = b0;
                    b0_s_old = b0_s;

                    double Upb   = (W.array() * resid_cur.array()     ).matrix().sum() / double(nobs)   + b0; //  / ;
                    double Upb_s = (W_s.array() * resid_s_cur.array() ).matrix().sum() / double(nobs_s) + b0_s; //  / ;

                    // double total_penalty = 0.0;
                    b0   = Upb; //soft_thresh(Upb, total_penalty);
                    b0_s = Upb_s; //soft_thresh(Upb_s, total_penalty);

                    //std::cout << "int U " << Upb << " int " << b0 << std::endl;

                    xbeta_cur.array() += (b0 - b0_old);
                    resid_cur.array() -= (b0 - b0_old);

                    xbeta_s_cur.array() += (b0_s - b0_s_old);
                    resid_s_cur.array() -= (b0_s - b0_s_old);
                }
                beta_old   = beta;
                beta_s_old = beta_s;

                for (int g = 0; g < ngroups; ++g)
                {
                    VectorXd beta_subs(2);
                    beta_subs(0) = beta(g);
                    beta_subs(1) = beta_s(g);

                    /*
                    VectorXd U_plus_beta = (x_list[g].transpose() *
                            ( W.array() * resid_cur.array()  ).matrix()).array() / double(nobs) +
                            eigenvals(g) * beta_subs.array();
                    */

                    VectorXd U_plus_beta(2);

                    U_plus_beta(0) = ( X.col(g).array() * W.array() * resid_cur.array()  ).matrix().sum() /
                        double(nobs) + eigenvals(g) * beta_subs(0);

                    U_plus_beta(1) = ( Xs.col(g).array() * W_s.array() * resid_s_cur.array()  ).matrix().sum() /
                        double(nobs_s) + eigenvals(g) * beta_subs(1);


                    U_plus_beta(0) *= mult_1;

                    double l1 = group_weights(g) * lam * tau;
                    double lgr = group_weights(g) * lam * (1.0 - tau);

                    // VectorXd beta_new = thresh_func(U_plus_beta, l1, gamma, l2, eigenvals(g));

                    VectorXd beta_new(2);
                    if (eigenvals(g) > 0.0)
                    {
                        //beta_new = thresh_func(U_plus_beta, penalty_adjustment, l1, gamma, l2, eigenvals(g));

                        if (tau > 0.0)
                        {
                            VectorXd beta_tmp(2);
                            for (int k = 0; k < 2; ++k)
                            {
                                double pencur = l1 * penalty_adjustment(k);
                                beta_tmp(k) = soft_thresh(U_plus_beta(k), pencur);
                            }

                            beta_new = thresh_func(beta_tmp, penalty_adjustment, lgr, gamma, lgr, eigenvals(g));
                        } else
                        {
                            beta_new = thresh_func(U_plus_beta, penalty_adjustment, lgr, gamma, lgr, eigenvals(g));
                        }

                        beta_new(0) *= mult_1;
                    } else
                    {
                        beta_new.setZero();
                    }

                    /*
                    if (g == 0)
                    {
                        std::cout << "U update " << U_plus_beta.transpose() << std::endl;
                        std::cout << "beta update " << beta_new.transpose() << std::endl;
                    }
                     */



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
                        //xbeta_cur += (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());
                        //resid_cur -= (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());

                        xbeta_cur.array() += (X.col(g).array() * (beta_new(0) - beta_subs(0)));
                        resid_cur.array() -= (X.col(g).array() * (beta_new(0) - beta_subs(0)));

                        xbeta_s_cur.array() += (Xs.col(g).array() * (beta_new(1) - beta_subs(1)));
                        resid_s_cur.array() -= (Xs.col(g).array() * (beta_new(1) - beta_subs(1)));
                    }

                }

                if (converged(beta, beta_old, tol) &&
                    converged(beta_s, beta_s_old, tol))
                {
                    break;
                }

            } // end BCD iter loop

            // if (converged_irls(deviance, deviance_old, tol_irls))
            // {
            //     niter(l) = ii + 1;
            //     break;
            // }
            if (converged(beta, beta_irls_old, tol_irls) &&
                converged(beta_s, beta_irls_s_old, tol_irls))
            {
                niter(l) = ii + 1;
                break;
            }
        } // end irls loop

        deviance_vec(l)   = deviance;
        deviance_s_vec(l) = deviance_s;

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
                        Named("niter")     = niter,
                        Named("lambda")    = lambda,
                        Named("tau")       = tau,
                        Named("deviance_z")  = deviance_vec,
                        Named("deviance_s")  = deviance_s_vec,
                        Named("penalty_adjustment") = penalty_adjustment,
                        Named("eigenvals") = eigenvals,
                        Named("family")    = family[0],
                        Named("penalty")   = penalty[0]);
}
