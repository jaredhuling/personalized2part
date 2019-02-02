
#include <Rcpp.h>
#include <RcppEigen.h>

#include "utils.h"

#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/DenseSymMatProd.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using namespace Spectra;


std::vector<MatrixXd > compute_xtx_list(std::vector<MatrixXd > & x_list,
                                        const Eigen::Map<Eigen::VectorXi> & groups,
                                        const int &ngroups,
                                        std::vector<std::vector<int> > & grp_idx,
                                        Eigen::VectorXd & weights)
{

    std::vector<MatrixXd > xtx_list(ngroups);

    int n = x_list[0].rows();

    for (int g = 0; g < ngroups; ++g)
    {
        std::vector<int> gr_idx = grp_idx[g];

        int numelem = gr_idx.size();

        MatrixXd XtWXtmp(MatrixXd(numelem, numelem).setZero().
                             selfadjointView<Eigen::Lower>().
                             rankUpdate( ((weights.array().sqrt().matrix()).asDiagonal() * x_list[g]).adjoint() ));

        XtWXtmp.array() /= n;

        xtx_list[g] = XtWXtmp;

    }

    return(xtx_list);
}

std::vector<MatrixXd > make_x_list(const Eigen::Map<Eigen::MatrixXd> & X,
                                   const Eigen::Map<Eigen::VectorXi> & groups,
                                   const int &ngroups,
                                   std::vector<std::vector<int> > & grp_idx)
{

    std::vector<MatrixXd > x_list(ngroups);

    int n = X.rows();


    for (int g = 0; g < ngroups; ++g)
    {
        std::vector<int> gr_idx = grp_idx[g];

        int numelem = gr_idx.size();

        MatrixXd sub(n, numelem);

        for (std::vector<int>::size_type v = 0; v < numelem; ++v)
        {
            int c_idx = gr_idx[v];
            sub.col(v) = X.col(c_idx);
        }

        x_list[g] = sub;

    }

    return(x_list);
}


Eigen::VectorXd compute_eigs(std::vector<MatrixXd > xtx_list)
{
    int ng = xtx_list.size();

    VectorXd eigvec(ng);

    for (int g = 0; g < ng; ++g)
    {
        DenseSymMatProd<double> op(xtx_list[g]);
        int ncv = 4;
        if (xtx_list[g].cols() < 4)
        {
            ncv = xtx_list[g].cols();
        }

        if (xtx_list[g].cols() > 1)
        {
            SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 1, ncv);

            eigs.init();
            eigs.compute(10000, 1e-10);
            VectorXd eigenvals = eigs.eigenvalues();

            eigvec(g) = eigenvals[0] * 1.0025; // multiply by an increasing factor to be safe
        } else
        {
            eigvec(g) = xtx_list[g](0,0);
        }
    }

    return(eigvec);
}

Eigen::VectorXd compute_eigs_twopart(const Eigen::Map<Eigen::MatrixXd> & X,
                                     const Eigen::Map<Eigen::MatrixXd> & Xs)
{
    int nc = X.cols();
    VectorXd eigvec(nc);

    int n  = X.rows();
    int ns = Xs.rows();

    for (int c = 0; c < nc; ++c)
    {
        double x_sq_1 = X.col(c).squaredNorm();
        double x_sq_2 = Xs.col(c).squaredNorm();

        eigvec(c) = std::max(x_sq_1 / double(n), x_sq_2 / double(ns));
        //eigvec(c) = std::max(x_sq_1, x_sq_2);
    }

    return(eigvec);
}

Eigen::VectorXd compute_eigs_twopart(const Eigen::Map<Eigen::MatrixXd> & X,
                                     const Eigen::Map<Eigen::MatrixXd> & Xs,
                                     Eigen::VectorXd & weights,
                                     Eigen::VectorXd & weights_s)
{
    int nc = X.cols();
    VectorXd eigvec(nc);

    int n  = X.rows();
    int ns = Xs.rows();

    VectorXd weights_sqrt = weights.array().sqrt();
    VectorXd weights_s_sqrt = weights_s.array().sqrt();

    for (int c = 0; c < nc; ++c)
    {
        double x_sq_1 = (X.col(c).array() * weights_sqrt.array()).matrix().squaredNorm();
        double x_sq_2 = (Xs.col(c).array() * weights_s_sqrt.array()).matrix().squaredNorm();

        eigvec(c) = std::max(x_sq_1 / double(n), x_sq_2 / double(ns));
    }

    return(eigvec);
}

Eigen::VectorXd setup_lambda(const Eigen::Map<Eigen::MatrixXd> & X,
                             std::vector<MatrixXd > & x_list,
                             const Eigen::Map<Eigen::VectorXd> & Y,
                             Eigen::VectorXd & weights,
                             const Eigen::Map<Eigen::VectorXi> & groups,
                             std::vector<std::vector<int> > & grp_idx,
                             Eigen::VectorXd & group_weights,
                             double b0,
                             U_func_ptr & U_func,
                             const int & nlambda,
                             const double & lambda_min_ratio,
                             std::string & penalty,
                             const double & alpha)
{
    //Eigen::VectorXd xty = X.transpose() * Y;
    int n = X.rows();
    int p = X.rows();

    VectorXd beta_init(p+1);
    beta_init(0) = b0;



    int ngroups = group_weights.size();

    VectorXd xbeta_cur = VectorXd::Constant(n, b0);

    VectorXd norms(ngroups);
    norms.setZero();

    for (int g = 0; g < ngroups; ++g)
    {

        // calculate U vector
        VectorXd U_vec = U_func(x_list[g], Y, weights, xbeta_cur, n).array();

        if (group_weights(g) > 0)
        {
            norms(g) = U_vec.norm() / group_weights(g);
        }

    }

    double lmax = norms.cwiseAbs().maxCoeff();






    //double lmax = xty.cwiseAbs().maxCoeff() / double(n);
    double lmin = lambda_min_ratio * lmax;

    VectorXd lambda_base(nlambda);

    lambda_base.setLinSpaced(nlambda, std::log(lmax), std::log(lmin));
    lambda_base = lambda_base.array().exp();

    return(lambda_base);
}



Eigen::VectorXd setup_lambda(const Eigen::Map<Eigen::MatrixXd> & X,
                             const Eigen::Map<Eigen::MatrixXd> & Xs,
                             const Eigen::Map<Eigen::VectorXd> & Z,
                             const Eigen::Map<Eigen::VectorXd> & S,
                             Eigen::VectorXd & weights,
                             Eigen::VectorXd & weights_s,
                             Eigen::VectorXd & group_weights,
                             double b0, double b0_s,
                             U_tp_func_ptr & U_func,
                             const int & nlambda,
                             const double & lambda_min_ratio,
                             std::string & penalty,
                             const double & alpha)
{
    //Eigen::VectorXd xty = X.transpose() * Y;
    int n = X.rows();
    int n_s = Xs.rows();
    int p = X.cols();

    VectorXd beta_init(p+1);
    beta_init(0) = b0;

    VectorXd beta_s_init(p+1);
    beta_s_init(0) = b0_s;



    int ngroups = group_weights.size();

    VectorXd xbeta_cur   = VectorXd::Constant(n, b0);
    VectorXd xbeta_s_cur = VectorXd::Constant(n_s, b0_s);

    VectorXd norms(ngroups);
    norms.setZero();

    VectorXd norms_z = norms;
    VectorXd norms_p = norms;

    MatrixXd U_mat(ngroups, 2);

    for (int g = 0; g < ngroups; ++g)
    {

        // calculate U vector
        VectorXd U_vec = U_func(X.col(g), Xs.col(g), Z, S, weights, weights_s,
                                xbeta_cur, xbeta_s_cur, n, n_s).array();

        norms_z(g) = std::abs(U_vec(0)) / group_weights(g);
        norms_p(g) = std::abs(U_vec(1)) / group_weights(g);

        U_mat.row(g) = U_vec;

    }


    double lmax_z = norms_z.cwiseAbs().maxCoeff();
    double lmax_p = norms_p.cwiseAbs().maxCoeff();

    VectorXd penalty_adjust(2);

    penalty_adjust(0) = lmax_z / (lmax_z + lmax_p);
    penalty_adjust(1) = lmax_p / (lmax_z + lmax_p);

    //penalty_adjust.array() /= penalty_adjust.cwiseAbs().minCoeff();

    for (int g = 0; g < ngroups; ++g)
    {
        if (group_weights(g) > 0)
        {
            norms(g) = U_mat.row(g).norm() / group_weights(g);
        }
    }


    double lmax = norms.cwiseAbs().maxCoeff();


    //double lmax = xty.cwiseAbs().maxCoeff() / double(n);
    double lmin = lambda_min_ratio * lmax;

    VectorXd lambda_base(nlambda);

    lambda_base.setLinSpaced(nlambda, std::log(lmax), std::log(lmin));
    lambda_base = lambda_base.array().exp();

    VectorXd lambda_plus_adjust(nlambda + 2);
    lambda_plus_adjust << penalty_adjust, lambda_base;

    return(lambda_plus_adjust);
}


Eigen::VectorXd setup_lambda(const Eigen::Map<Eigen::MatrixXd> & X,
                             const Eigen::Map<Eigen::VectorXd> & Y,
                             Eigen::VectorXd & weights,
                             const Eigen::Map<Eigen::VectorXi> & groups,
                             std::vector<std::vector<int> > & grp_idx,
                             Eigen::VectorXd & group_weights,
                             const int & nlambda,
                             const double & lambda_min_ratio,
                             std::string & penalty,
                             const double & alpha)
{
    Eigen::VectorXd xty = X.transpose() * (weights.array() * Y.array()).matrix();
    int n = X.rows();
    //int p = X.cols();

    int ngroups = group_weights.size();

    VectorXd norms(ngroups);
    norms.setZero();


    for (int g = 0; g < ngroups; ++g)
    {
        if (group_weights(g) > 0)
        {
            std::vector<int> gr_idx = grp_idx[g];

            int numelem = gr_idx.size();

            double norm_cur = 0.0;

            for (std::vector<int>::size_type v = 0; v < numelem; ++v)
            {
                int c_idx = gr_idx[v];
                norm_cur += std::pow(xty(c_idx), 2);
            }

            norms(g) = std::sqrt(norm_cur) / group_weights(g);
        }
    }

    double lmax = norms.cwiseAbs().maxCoeff() / double(n);
    double lmin = lambda_min_ratio * lmax;

    VectorXd lambda_base(nlambda);

    lambda_base.setLinSpaced(nlambda, std::log(lmax), std::log(lmin));
    lambda_base = lambda_base.array().exp();

    return(lambda_base);
}

std::vector<std::vector<int> > get_group_indexes(const VectorXi &groups,
                                                 const VectorXi &unique_groups,
                                                 const int &ngroups,
                                                 const int &nvars)
{
    std::vector<std::vector<int> > grp_idx(ngroups);
    //grp_idx.reserve(ngroups);

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
    // if group weights were not specified,
    // then set the group weight for each
    // group to be the sqrt of the size of the
    // group
    /*
    if (default_group_weights)
    {
        group_weights.resize(ngroups);
        for (int g = 0; g < ngroups; ++g)
        {
            group_weights(g) = std::sqrt(double(grp_idx[g].size()));
        }
    }
     */
    return(grp_idx);
}


bool converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance)
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

bool converged_irls(double deviance, double deviance_prev, const double& tolerance)
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
