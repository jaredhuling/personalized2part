
#include "fit_gamma.h"


double gammamod::grad_func(const VectorXd &x_col,
                             VectorXd &xbeta)
{

    double U;

    // gamma part
    U = (x_col.array() *
        (weights.array() * (Y.array() * (-1.0 * xbeta.array()).exp() - 1.0  ))).matrix().sum() / double(nobs);

    return(U);
}


double gammamod::grad_func(int col_idx)
{
    double U;

    // gamma part
    U = (X.col(col_idx).array() *
        (weights.array() * (Y.array() * (-1.0 * xbeta_cur.array()).exp() - 1.0  ))).matrix().sum() / double(nobs);

    return(U);
}

void gammamod::update_strongrule(int lam_idx)
{
    double lam_prev = 0.0;
    double lam_cur = lambda(lam_idx);
    if (lam_idx > 0)
    {
        lam_prev = lambda(lam_idx - 1);
    }

    active_set.setZero();

    double grad_cur;
    double threshed_grad_cur;

    for (int j = 0; j < nvars; ++j)
    {
        if (group_weights(j) > 0.0)
        {
            grad_cur = grad_func(j);

            double grp_thresh = (1.0 - tau) * group_weights(j) * (2.0 * lam_cur - lam_prev);
            double individ_thresh = tau * group_weights(j) * (2.0 * lam_cur - lam_prev);

            double uval = grad_cur;
            threshed_grad_cur = soft_thresh(uval, individ_thresh);

            if (std::abs(threshed_grad_cur) >= grp_thresh)
            {
                //std::cout << "grp norm: " << grad_cur.norm() << " thresh: " << grp_thresh << " actset: " << active_set(j) << std::endl;
                active_set(j) = 1;
            }

        } else
        {
            active_set(j) = 1;
        }
    }

    // std::cout << "lam " << lam_idx << " num active init: " << active_set.sum() << std::endl;
}

void gammamod::check_kkt(int lam_idx)
{
    any_violations = false;
    double lam = lambda(lam_idx);
    double grad_cur;
    double threshed_grad_cur;

    int num_violations = 0;

    for (int j = 0; j < nvars; ++j)
    {
        if (active_set(j) == 0)
        {
            if (group_weights(j) > 0.0)
            {
                double l1  = group_weights(j) * lam * tau;
                double lgr = group_weights(j) * lam * (1.0 - tau);

                grad_cur = grad_func(j);
                double uval = grad_cur;
                threshed_grad_cur = soft_thresh(uval, l1);

                if (std::abs(threshed_grad_cur) >= lgr)
                {
                    any_violations = true;
                    active_set(j) = 1;
                    ++num_violations;
                }

            }
        }
    }

    // std::cout << "lam " << lam_idx << " num violations: " << num_violations << " num active: " << active_set.sum() << std::endl;
}



// phi_j(v) function for cooperative lasso
VectorXd gammamod::phi_j_v(VectorXd & v, int & j)
{
    int vlen = v.size();
    VectorXd retvec(vlen);

    retvec.setZero();

    if (v(j) > 0.0)
    {
        // calculate v plus
        for (int k = 0; k < vlen; ++k)
        {
            retvec(k) = std::max(0.0, v(k));
        }
    } else if (v(j) < 0.0)
    {
        // calculate v minus
        for (int k = 0; k < vlen; ++k)
        {
            retvec(k) = std::max(0.0, -1.0 * v(k));
        }
    }

    return(retvec);
}


bool gammamod::converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance)
{
    for (int i = 0; i < cur.rows(); i++)
    {
        if ( (std::abs(cur(i)) > 1e-13 && std::abs(prev(i)) <= 1e-13) ||
             (std::abs(cur(i)) <= 1e-13 && std::abs(prev(i)) > 1e-13) ) {
            return 0;
        }
        if (std::abs(cur(i)) > 1e-13 && std::abs(prev(i)) > 1e-13 &&
            std::pow( (cur(i) - prev(i)) / prev(i), 2) > tolerance) {
            return 0;
        }
    }
    return 1;
}

bool gammamod::converged_irls(double deviance, double deviance_prev, const double& tolerance)
{
    //return (stopRule(beta, beta_prev_irls, tol_irls));

    if (std::abs(deviance - deviance_prev) / (0.1 + std::abs(deviance)) < tolerance)
    {
        return true;
    } else
    {
        return false;
    }

}

VectorXd gammamod::compute_eigs_twopart()
{
    VectorXd eigvec(nvars);

    VectorXd weights_sqrt   = W.array().sqrt();

    for (int c = 0; c < nvars; ++c)
    {
        double x_sq_1 = (X.col(c).array() * weights_sqrt.array()).matrix().squaredNorm();

        eigvec(c) = x_sq_1 / double(nobs);
    }


    return(eigvec);
}

double gammamod::soft_thresh(double & a, double & lambda)
{
    //double thresh_fact = std::max(0.0, 1.0 - lambda / std::abs(a) );
    //double retval = a * thresh_fact;

    double retval = a * std::max(0.0, 1.0 - lambda / std::abs(a) );

    return(retval);
}

VectorXd gammamod::block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    VectorXd beta_tmp(alen);
    if (l1 > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            double l1_apply = penalty_factor(j) * l1;
            beta_tmp(j) = soft_thresh(a(j), l1_apply);
        }
    } else
    {
        beta_tmp = a;
    }

    double anorm = beta_tmp.norm();
    double thresh_fact;

    if (denom > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            //thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / anorm) / denom;
            thresh_fact = std::max(0.0, 1.0 - lambda / anorm) / denom;
            retval(j) = beta_tmp(j) * thresh_fact;
        }
    } else
    {
        retval.setZero();
    }

    return(retval);
}

/*
VectorXd gammamod::coop_block_soft_thresh_tp(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    VectorXd beta_tmp(alen);
    if (l1 > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            double l1_apply = penalty_factor(j) * l1;
            beta_tmp(j) = soft_thresh(a(j), l1_apply);
        }
    } else
    {
        beta_tmp = a;
    }

    double thresh_fact;

    if (denom > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            VectorXd phi_j_vec = phi_j_v(beta_tmp, j);
            //thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / phi_j_vec.norm());
            thresh_fact = std::max(0.0, 1.0 - lambda / phi_j_vec.norm());
            retval(j) = beta_tmp(j) * thresh_fact / denom;
        }
    } else
    {
        retval.setZero();
    }

    return(retval);
}

*/

/*
VectorXd twopart::thresh_func(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
{
    if (penalty == "grp.lasso")
    {
        return(coop_block_soft_thresh_tp(a, penalty_factor, lambda, l1, denom));
    } else
    {
        return(block_soft_thresh_tp(a, penalty_factor, lambda, l1, denom));
    }
}*/

void gammamod::initialize()
{

    b0   = 0.0;

    if (intercept)
    {
        double sbar = (Y.array() * weights.array()).matrix().sum() / weights.sum(); // double(nobs_s);
        //double sbar = S.sum() / double(nobs_s);
        b0 = log(sbar);
    }

    //b0 = Y.sum() / double(nobs);
    b0_old   = b0;

    // ----------------------------------------- //
    // ------        Set up groups        ------ //
    // ----------------------------------------- //

    set_up_groups();


    // ----------------------------------------- //
    // ------        Set up lambda        ------ //
    // ----------------------------------------- //

    // run once to set scale of likelihoods, again to double check
    set_up_lambda();


    set_up_lambda();



    if (intercept)
    {
        xbeta_cur.array() += b0;
    }


    // ----------------------------------------- //
    // ---    Set up threshold function   ------ //
    // ----------------------------------------- //


    thresh_func = &gammamod::block_soft_thresh;

    /*
    if (penalty == "grp.lasso")
    {
        thresh_func = &gammamod::block_soft_thresh;
    } else if (penalty == "coop.lasso")
    {
        thresh_func = &gammamod::coop_block_soft_thresh_tp;
    } else
    {
        thresh_func = &gammamod::block_soft_thresh;
    }*/


}

void gammamod::set_up_lambda()
{
    if (need_to_compute_lambda)
    {
        VectorXd xbeta_init   = VectorXd::Constant(nobs, b0);

        VectorXd norms(ngroups);
        norms.setZero();

        VectorXd U_mat(ngroups);

        for (int g = 0; g < ngroups; ++g)
        {

            // calculate U vector
            double U = grad_func(X.col(g), xbeta_init);

            if (group_weights(g) > 0.0)
            {
                norms(g) = std::abs(U) / group_weights(g);
            }

            U_mat(g) = U;

        }


        double lmax = norms.cwiseAbs().maxCoeff();


        //penalty_adjust.array() /= penalty_adjust.cwiseAbs().minCoeff();

        if (penalty == "grp.lasso")
        {
            for (int g = 0; g < ngroups; ++g)
            {
                if (group_weights(g) > 0.0)
                {
                    norms(g) = std::abs(U_mat(g)) / group_weights(g);
                }
            }
            lmax = norms.cwiseAbs().maxCoeff();
        }

        VectorXd penalty_adjust(2);



        //std::cout << lmax_z  << " " << lmax_p << std::endl;
        //std::cout << penalty_adjust.transpose() << std::endl;

        //double lmax = xty.cwiseAbs().maxCoeff() / double(n);
        double lmin = lambda_min_ratio * lmax;

        VectorXd lambda_base(nlambda);

        lambda_base.setLinSpaced(nlambda, std::log(lmax), std::log(lmin));
        lambda_base = lambda_base.array().exp();

        lambda = lambda_base;

        //penalty_adjustment.array() = 1.0;

        //penalty_adjustment(0) = 0.25;
    } else
    {
        lambda = lambda_given;
    }


}

void gammamod::set_up_groups()
{
    bool default_group_weights = bool(group_weights_given.size() < 1);

    // set up groups
    //std::vector<std::vector<int> > grp_idx = get_group_indexes(groups, unique_groups, ngroups, nvars);

    //std::vector<std::vector<int> > grp_idx(ngroups);

    for (int g = 0; g < ngroups; ++g)
    {
        // find all variables in group number g
        std::vector<int> idx_tmp;
        for (int v = 0; v < nvars; ++v)
        {
            if (groups(v) == unique_groups(g))
            {
                idx_tmp.push_back(v);
            }
        }

        grp_idx[g] = idx_tmp;
    }

    if (default_group_weights)
    {
        group_weights.resize(ngroups);
        for (int g = 0; g < ngroups; ++g)
        {
            group_weights(g) = std::sqrt(double(grp_idx[g].size()));
        }
    } else
    {
        group_weights = group_weights_given;
    }
}

VectorXi gammamod::fit_path()
{
    VectorXi itervec(nlambda);

    itervec.array() = maxit_irls;

    // std::cout << itervec.size() << deviance_vec.size() << deviance_s_vec.size() << std::endl;
    if (!strongrule)
    {
        active_set.array() = 1;
    }

    // loop over lamba values
    for (int l = 0; l < nlambda; ++l)
    {
        double lam = lambda(l);
        deviance   = 1e30;

        if (strongrule)
        {
            update_strongrule(l);
        }

        for (int ct = 0; ct < 5; ++ct)
        {
            // irls iters
            for (int ii = 0; ii < maxit_irls; ++ii)
            {
                beta_irls_old   = beta;
                deviance_old    = deviance;

                // -------------------------------------------------- //
                //                update quad approx                  //


                /////// ----------------
                ///////   GAMMA UPDATE
                /////// ----------------

                mu = Y.array() * (-1.0 * xbeta_cur.array()).exp();
                W = weights.array() * mu.array();

                resid_cur = (mu.array() - 1.0) / mu.array();

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

                deviance = (Y.array().log() - xbeta_cur.array() - mu.array() + 1.0).matrix().sum();
                deviance *= -2.0;

                if (converged_irls(deviance,   deviance_old,   tol_irls))
                {
                    itervec(l) = ii;
                    break;
                }

                //Xsq = (W.array().sqrt().matrix().asDiagonal() * datX).array().square().colwise().sum();

                //resid_cur = Y - xbeta_cur;


                //VectorXd realweights = W.array() * weights.array();
                //xtx_list = compute_xtx_list(x_list, groups, ngroups, grp_idx, W);

                // compute largest eigenvalues within each group
                //eigenvals = compute_eigs(xtx_list);
                VectorXd eigenvals = compute_eigs_twopart();

                //std::cout << "eigs " << eigenvals.head(5).transpose() << std::endl;

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
                        b0 = Upb; //soft_thresh(Upb_s, total_penalty);

                        xbeta_cur.array() += (b0 - b0_old);
                        resid_cur.array() -= (b0 - b0_old);
                    }
                    beta_old   = beta;

                    for (int g = 0; g < ngroups; ++g)
                    {
                        if (active_set(g) == 1)
                        {
                            double beta_subs = beta(g);

                            /*
                            VectorXd U_plus_beta = (x_list[g].transpose() *
                            ( W.array() * resid_cur.array()  ).matrix()).array() / double(nobs) +
                            eigenvals(g) * beta_subs.array();
                            */

                            VectorXd U_plus_beta(1);

                            U_plus_beta(0) = ( X.col(g).array() * W.array() * resid_cur.array()  ).matrix().sum() /
                                double(nobs) + eigenvals(g) * beta_subs;

                            //U_plus_beta(0) *= mult_1;

                            double l1  = group_weights(g) * lam * tau;
                            double lgr = group_weights(g) * lam * (1.0 - tau);

                            // this is how we call a pointer to a member function from within the class
                            VectorXd beta_new = (this->*thresh_func)(U_plus_beta, penalty_adjustment, lgr, l1, eigenvals(g));
                            //VectorXd beta_new = thresh_func(U_plus_beta, penalty_adjustment, lgr, l1, eigenvals(g));

                            /*
                            if (g == 0)
                            {
                            std::cout << "U update " << U_plus_beta.transpose() << std::endl;
                            std::cout << "beta update " << beta_new.transpose() << std::endl;
                            }
                            */



                            bool anychanged = false;
                            if (beta_subs != beta_new(0))
                            {
                                anychanged = true;
                                beta(g)    = beta_new(0);
                            }
                            /*
                            for (int k = 0; k < 2; ++k)
                            {
                                if (beta_subs(k) != beta_new(k))
                                {
                                    anychanged = true;
                                    beta(g)    = beta_new(0);
                                    beta_s(g)  = beta_new(1);
                                }
                            }*/


                            // update residual if any estimate changed
                            if (anychanged)
                            {
                                //xbeta_cur += (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());
                                //resid_cur -= (x_list[g] * (beta_new.array() - beta_subs.array()).matrix());

                                xbeta_cur.array() += (X.col(g).array() * (beta_new(0) - beta_subs));
                                resid_cur.array() -= (X.col(g).array() * (beta_new(0) - beta_subs));
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
                if (converged(beta, beta_irls_old, tol_irls))
                {
                    //std::cout << "lam idx: " << l << std::endl;
                    itervec(l) = ii + 1;
                    break;
                }
            } // end irls loop

            if (strongrule)
            {
                check_kkt(l);
                if (!any_violations)
                {
                    break;
                }
            } else
            {
                break;
            }

        }

        deviance_vec(l)   = deviance;

        beta_mat.col(l).tail(nvars)   = beta;
        if (intercept)
        {
            beta_mat(0, l)   = b0;
        }
    }

    return(itervec);
}
