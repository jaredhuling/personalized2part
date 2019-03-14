#ifndef FITGAMMA_H
#define FITGAMMA_H
#include "RcppEigen.h"
#include "params.h"


using namespace Rcpp;
using namespace RcppEigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Map;

typedef Map<MatrixXd> MapMatd;
typedef Map<const MatrixXd> cMapMatd;
typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;
typedef Map<const VectorXd> cMapVecd;
typedef Map<const VectorXi> cMapVeci;


class gammamod
{
    protected:


        // pointer we will set to one of the thresholding functions
        typedef VectorXd (gammamod::*thresh_func_ptr)(VectorXd &value, VectorXd & penalty_factor,
                          double &penalty, double &l1, double &denom);



        const cMapMatd X;
        const cMapVecd Y;
        const cMapVecd lambda_given;
        const cMapVecd weights;
        const cMapVecd offset;
        const cMapVeci groups, unique_groups;
        const cMapVecd group_weights_given;
        params P;

        double tau;
        int maxit, maxit_irls;
        double tol, tol_irls;

        bool intercept;
        std::string penalty;
        bool strongrule;

        bool need_to_compute_lambda;
        int nlambda;
        double lambda_min_ratio;

        int ngroups;
        std::vector<std::vector<int> > grp_idx;
        VectorXd eigenvals;


        int nobs, nvars;

        VectorXd beta, beta_old, beta_irls_old;
        MatrixXd beta_mat;

        VectorXd mu, xbeta_cur, resid_cur, W;
        VectorXd deviance_vec;
        VectorXi active_set;

        bool any_violations = false;

        double b0, b0_old, deviance, deviance_old;

        bool scale_set = false;

        VectorXd penalty_adjustment, ZZ, lambda, group_weights;

        thresh_func_ptr thresh_func;
        //VectorXd thresh_func(VectorXd &value, VectorXd & penalty_factor, double &penalty, double &l1, double &denom);


        double grad_func(const VectorXd &x_col,
                         VectorXd &xbeta);

        virtual void set_up_lambda();
        virtual void set_up_groups();
        virtual double grad_func(int col_idx);
        virtual void update_strongrule(int lam_idx);
        virtual void check_kkt(int lam_idx);

        // phi_j(v) function for cooperative lasso
        VectorXd phi_j_v(VectorXd & v, int & j);

        bool converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance);

        bool converged_irls(double deviance, double deviance_prev, const double& tolerance);

        VectorXd compute_eigs_twopart();

        double soft_thresh(double & a, double & lambda);

        VectorXd block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom);
        //VectorXd coop_block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom);




    public:
        gammamod(const Eigen::Ref<const MatrixXd> & X_,
                 const Eigen::Ref<const VectorXd> & Y_,
                 const Eigen::Ref<const VectorXd> & weights_,
                 const Eigen::Ref<const VectorXd> & offset_,
                 const Eigen::Ref<const VectorXd> & lambda_given_,
                 const Eigen::Ref<const VectorXi> & groups_,
                 const Eigen::Ref<const VectorXi> & unique_groups_,
                 const Eigen::Ref<const VectorXd> & group_weights_given_,
                 params & P_) :
        X(X_.data(), X_.rows(), X_.cols()),
        Y(Y_.data(), Y_.size()),
        lambda_given(lambda_given_.data(), lambda_given_.size()),
        weights(weights_.data(), weights_.size()),
        offset(offset_.data(), offset_.size()),
        groups(groups_.data(), groups_.size()),
        unique_groups(unique_groups_.data(), unique_groups_.size()),
        group_weights_given(group_weights_given_.data(), group_weights_given_.size()),
        P(P_),
        tau(P.tau),
        maxit(P.maxit),
        maxit_irls(P.maxit_irls),
        tol(P.tol),
        tol_irls(P.tol_irls),
        intercept(P.intercept_s),
        penalty(P.penalty),
        strongrule(P.strongrule),
        need_to_compute_lambda(lambda_given.size() < 1),
        nlambda((need_to_compute_lambda) ? P.nlambda : lambda_given.size()),
        lambda_min_ratio(P.lambda_min_ratio),
        ngroups(unique_groups.size()),
        grp_idx(ngroups),
        eigenvals(ngroups),
        nobs(X.rows()),
        nvars(X.cols()),
        beta(VectorXd::Zero(nvars)),
        beta_old(VectorXd::Zero(nvars)),
        beta_irls_old(VectorXd::Zero(nvars)),
        beta_mat(MatrixXd::Zero(nvars+1, nlambda)),
        mu(nobs), xbeta_cur(offset), resid_cur(nobs),
        W(nobs),
        deviance_vec(nlambda),
        active_set(VectorXi::Zero(nvars))
        {}


        virtual void initialize();
        virtual VectorXi fit_path();

        virtual MatrixXd get_beta()    { return beta_mat; }
        virtual VectorXd get_lambda()  { return lambda; }
        virtual VectorXd get_dev()     { return deviance_vec; }
        virtual VectorXd get_pen_adj() { return penalty_adjustment; }
        virtual VectorXd get_eigs()    { return eigenvals; }


};

#endif
