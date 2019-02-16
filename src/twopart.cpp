
#include "twopart.h"


VectorXd twopart::grad_func(const VectorXd &x_col,
                            const VectorXd &x_col_s,
                            VectorXd &xbeta,
                            VectorXd &xbeta_s)
{
    VectorXd U_vec(2);

    // binomial part
    U_vec(0) = mult_1 * (x_col.array() *
        (weights.array() * (Z.array() / (1.0 + (Z.array() * xbeta.array()).exp()  )))).matrix().sum() / double(nobs);

    // gamma part
    U_vec(1) = (x_col_s.array() *
        (weights_s.array() * (S.array() * (-1.0 * xbeta_s.array()).exp() - 1.0  ))).matrix().sum() / double(nobs_s);

    return(U_vec);
}


VectorXd twopart::grad_func(int col_idx)
{
    VectorXd U_vec(2);

    // binomial part
    U_vec(0) = mult_1 * (X.col(col_idx).array() *
        (weights.array() * (Z.array() / (1.0 + (Z.array() * xbeta_cur.array()).exp()  )))).matrix().sum() / double(nobs);

    // gamma part
    U_vec(1) = (Xs.col(col_idx).array() *
        (weights_s.array() * (S.array() * (-1.0 * xbeta_s_cur.array()).exp() - 1.0  ))).matrix().sum() / double(nobs_s);

    return(U_vec);
}

void twopart::update_strongrule(int lam_idx)
{
    double lam_prev = 0.0;
    double lam_cur = lambda(lam_idx);
    if (lam_idx > 0)
    {
        lam_prev = lambda(lam_idx - 1);
    }

    active_set.setZero();

    VectorXd grad_cur(2);
    VectorXd threshed_grad_cur(2);

    for (int j = 0; j < nvars; ++j)
    {
        if (group_weights(j) > 0.0)
        {
            grad_cur = grad_func(j);

            double grp_thresh = (1.0 - tau) * group_weights(j) * (2.0 * lam_cur - lam_prev);
            double individ_thresh = tau * group_weights(j) * (2.0 * lam_cur - lam_prev);

            for (int g = 0; g < 2; ++g)
            {
                threshed_grad_cur(g) = soft_thresh(grad_cur(g), individ_thresh);
            }

            if (penalty == "grp.lasso")
            {
                if (threshed_grad_cur.norm() >= grp_thresh)
                {
                    //std::cout << "grp norm: " << grad_cur.norm() << " thresh: " << grp_thresh << " actset: " << active_set(j) << std::endl;
                    active_set(j) = 1;
                }
                /*else
                {
                    if (tau > 0.0)
                    {
                        for (int g = 0; g < grad_cur.size(); ++g)
                        {
                            if (std::abs(grad_cur(g)) >= individ_thresh)
                            {
                                active_set(j) = 1;
                                break;
                            }
                        }
                    }
                }*/
            } else
            {
                for (int g = 0; g < grad_cur.size(); ++g)
                {
                    /*
                    if (tau > 0.0)
                    {
                        if (std::abs(grad_cur(g)) >= individ_thresh)
                        {
                            active_set(j) = 1;
                            break;
                        }
                    }*/
                    VectorXd phi_j_vec = phi_j_v(threshed_grad_cur, g);
                    if (phi_j_vec.norm() >= grp_thresh)
                    {
                        active_set(j) = 1;
                        break;
                    }
                }
            }
        } else
        {
            active_set(j) = 1;
        }
    }

    // std::cout << "lam " << lam_idx << " num active init: " << active_set.sum() << std::endl;
}

void twopart::check_kkt(int lam_idx)
{
    any_violations = false;
    double lam = lambda(lam_idx);
    VectorXd grad_cur(2);
    VectorXd threshed_grad_cur(2);

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
                for (int g = 0; g < 2; ++g)
                {
                    threshed_grad_cur(g) = soft_thresh(grad_cur(g), l1);
                }

                if (penalty == "grp.lasso")
                {
                    if (threshed_grad_cur.norm() >= lgr)
                    {
                        any_violations = true;
                        active_set(j) = 1;
                        ++num_violations;
                    }
                } else
                {
                    for (int g = 0; g < 2; ++g)
                    {
                        VectorXd phi_j_vec = phi_j_v(threshed_grad_cur, g);
                        if (phi_j_vec.norm() >= lgr)
                        {
                            any_violations = true;
                            active_set(j) = 1;
                            ++num_violations;
                            break;
                        }
                    }
                }
            }
        }
    }

    // std::cout << "lam " << lam_idx << " num violations: " << num_violations << " num active: " << active_set.sum() << std::endl;
}



// phi_j(v) function for cooperative lasso
VectorXd twopart::phi_j_v(VectorXd & v, int & j)
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


bool twopart::converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance)
{
    for (unsigned i = 0; i < cur.rows(); i++)
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

bool twopart::converged_irls(double deviance, double deviance_prev, const double& tolerance)
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

VectorXd twopart::compute_eigs_twopart()
{
    VectorXd eigvec(nvars);

    VectorXd weights_sqrt = W.array().sqrt();
    VectorXd weights_s_sqrt = W_s.array().sqrt();

    for (int c = 0; c < nvars; ++c)
    {
        double x_sq_1 = (X.col(c).array() * weights_sqrt.array()).matrix().squaredNorm();
        double x_sq_2 = (Xs.col(c).array() * weights_s_sqrt.array()).matrix().squaredNorm();

        eigvec(c) = std::max(x_sq_1 / double(nobs), x_sq_2 / double(nobs_s));
    }


    return(eigvec);
}

double twopart::soft_thresh(double & a, double & lambda)
{
    //double thresh_fact = std::max(0.0, 1.0 - lambda / std::abs(a) );
    //double retval = a * thresh_fact;

    double retval = a * std::max(0.0, 1.0 - lambda / std::abs(a) );

    return(retval);
}

VectorXd twopart::block_soft_thresh_tp(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    VectorXd beta_tmp(alen);
    if (l1 > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            beta_tmp(j) = soft_thresh(a(j), l1);
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
            thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / anorm) / denom;
            retval(j) = beta_tmp(j) * thresh_fact;
        }
    } else
    {
        retval.setZero();
    }

    return(retval);
}

VectorXd twopart::coop_block_soft_thresh_tp(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    VectorXd beta_tmp(alen);
    if (l1 > 0.0)
    {
        for (int j = 0; j < alen; ++j)
        {
            beta_tmp(j) = soft_thresh(a(j), l1);
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
            thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / phi_j_vec.norm());
            retval(j) = beta_tmp(j) * thresh_fact / denom;
        }
    } else
    {
        retval.setZero();
    }

    return(retval);
}

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

void twopart::initialize()
{
    mult_1 = 1.0;
    if (opposite_signs)
    {
        mult_1 = -1.0;
    }

    ZZ.resize(nobs);
    ZZ = (Z.array() + 1) / 2;

    b0   = 0.0;
    b0_s = 0.0;
    if (intercept)
    {
        double zbar = ( (ZZ.array() * weights.array()) ).matrix().sum() /  weights.sum(); // double(nobs);
        b0 = log(zbar / (1.0 - zbar));

        double sbar = (S.array() * weights_s.array()).matrix().sum() / weights_s.sum(); // double(nobs_s);
        //double sbar = S.sum() / double(nobs_s);
        b0_s = log(sbar);
    }

    //b0 = Y.sum() / double(nobs);
    b0_old   = b0;
    b0_s_old = b0_s;



    // ----------------------------------------- //
    // ------        Set up groups        ------ //
    // ----------------------------------------- //

    set_up_groups();


    // ----------------------------------------- //
    // ------        Set up lambda        ------ //
    // ----------------------------------------- //

    set_up_lambda();


    // set up linear predictors
    if (intercept)
    {
        xbeta_cur.array()   += b0;
        xbeta_s_cur.array() += b0_s;
    } else
    {
        //xbeta_cur.setZero();
        //xbeta_s_cur.setZero();
    }


    // ----------------------------------------- //
    // ---    Set up threshold function   ------ //
    // ----------------------------------------- //


    if (penalty == "grp.lasso")
    {
        thresh_func = &twopart::block_soft_thresh_tp;
    } else if (penalty == "coop.lasso")
    {
        thresh_func = &twopart::coop_block_soft_thresh_tp;
    } else
    {
        thresh_func = &twopart::block_soft_thresh_tp;
    }


}

void twopart::set_up_lambda()
{
    if (need_to_compute_lambda)
    {
        VectorXd xbeta_init   = VectorXd::Constant(nobs, b0);
        VectorXd xbeta_s_init = VectorXd::Constant(nobs_s, b0_s);

        VectorXd norms(ngroups);
        norms.setZero();

        VectorXd norms_z = norms;
        VectorXd norms_p = norms;

        MatrixXd U_mat(ngroups, 2);

        for (int g = 0; g < ngroups; ++g)
        {

            // calculate U vector
            VectorXd U_vec = grad_func(X.col(g), Xs.col(g), xbeta_init, xbeta_s_init);

            if (group_weights(g) > 0.0)
            {
                if (penalty == "grp.lasso")
                {
                    norms_z(g) = std::abs(U_vec(0)) / group_weights(g);
                    norms_p(g) = std::abs(U_vec(1)) / group_weights(g);
                } else if (penalty == "coop.lasso")
                {
                    VectorXd norms_phis(2);
                    VectorXd phi_j_vec(2);

                    for (int j = 0; j < 2; ++j)
                    {
                        phi_j_vec = phi_j_v(U_vec, j);
                        norms_phis(j) = phi_j_vec.norm() / group_weights(g);
                    }

                    norms_z(g) = norms_phis(0);
                    norms_p(g) = norms_phis(1);
                }
            }

            U_mat.row(g) = U_vec;

        }


        double lmax_z = norms_z.cwiseAbs().maxCoeff();
        double lmax_p = norms_p.cwiseAbs().maxCoeff();


        //penalty_adjust.array() /= penalty_adjust.cwiseAbs().minCoeff();

        for (int g = 0; g < ngroups; ++g)
        {
            if (group_weights(g) > 0.0)
            {
                norms(g) = U_mat.row(g).norm() / group_weights(g);
            }
        }


        double lmax = norms.cwiseAbs().maxCoeff();

        VectorXd penalty_adjust(2);

        penalty_adjust(0) = lmax_z / (lmax_z + lmax_p);
        penalty_adjust(1) = lmax_p / (lmax_z + lmax_p);


        //double lmax = xty.cwiseAbs().maxCoeff() / double(n);
        double lmin = lambda_min_ratio * lmax;

        VectorXd lambda_base(nlambda);

        lambda_base.setLinSpaced(nlambda, std::log(lmax), std::log(lmin));
        lambda_base = lambda_base.array().exp();

        lambda = lambda_base;

        penalty_adjustment = penalty_adjust;
        penalty_adjustment.array() /= penalty_adjustment.maxCoeff();

        penalty_adjustment.array() = 1.0;
    } else
    {
        lambda = lambda_given;
        penalty_adjustment.resize(2);
        penalty_adjustment.array() = 1.0;
    }


}

void twopart::set_up_groups()
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

VectorXi twopart::fit_path()
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
        deviance_s = 1e30;

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

                beta_irls_s_old = beta_s;
                deviance_s_old  = deviance_s;

                // -------------------------------------------------- //
                //                update quad approx                  //


                // calculate mean function
                mu = 1.0 / (1.0 + (-1.0 * xbeta_cur.array()).exp() );

                //std::cout << "min prob" << p.minCoeff() << std::endl;


                // construct weights and multiply by user-specified weights
                //W = weights.array() * p.array() * (1.0 - p.array());
                VectorXd var_vec = mu.array() * (1.0 - mu.array());
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
                        if (mu(ii) > 1e-5)
                        {
                            deviance -= std::log(mu(ii));
                        } else
                        {
                            // don't divide by zero
                            deviance -= std::log(1e-5);
                        }

                    } else
                    {
                        if (mu(ii) <= 1.0 - 1e-5)
                        {
                            deviance -= std::log((1.0 - mu(ii)));
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
                resid_cur = (ZZ.array() - mu.array()) / var_vec.array(); // + xbeta_cur.array() * W.array().sqrt();

                /////// ----------------
                ///////   GAMMA UPDATE
                /////// ----------------

                mu_s = S.array() * (-1.0 * xbeta_s_cur.array()).exp();
                W_s = weights_s.array() * mu_s.array();

                resid_s_cur = (mu_s.array() - 1.0) / mu_s.array();

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

                deviance_s = (S.array().log() - xbeta_s_cur.array() - mu_s.array() + 1.0).matrix().sum();
                deviance_s *= -2.0;

                if (converged_irls(deviance,   deviance_old,   tol_irls) &&
                    converged_irls(deviance_s, deviance_s_old, tol_irls))
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
                        if (active_set(g) == 1)
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

                            U_plus_beta(0) = mult_1 * ( X.col(g).array() * W.array() * resid_cur.array()  ).matrix().sum() /
                                double(nobs) + eigenvals(g) * beta_subs(0);

                            U_plus_beta(1) = ( Xs.col(g).array() * W_s.array() * resid_s_cur.array()  ).matrix().sum() /
                                double(nobs_s) + eigenvals(g) * beta_subs(1);


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

                                xbeta_cur.array() += mult_1 * (X.col(g).array() * (beta_new(0) - beta_subs(0)));
                                resid_cur.array() -= mult_1 * (X.col(g).array() * (beta_new(0) - beta_subs(0)));

                                xbeta_s_cur.array() += (Xs.col(g).array() * (beta_new(1) - beta_subs(1)));
                                resid_s_cur.array() -= (Xs.col(g).array() * (beta_new(1) - beta_subs(1)));
                            }
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
        deviance_s_vec(l) = deviance_s;

        beta_mat.col(l).tail(nvars)   = beta;
        beta_s_mat.col(l).tail(nvars) = beta_s;
        if (intercept)
        {
            beta_mat(0, l)   = b0;
            beta_s_mat(0, l) = b0_s;
        }
    }

    return(itervec);
}
