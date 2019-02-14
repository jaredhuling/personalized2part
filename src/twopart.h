#ifndef TWOPART_H
#define TWOPART_H
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


class twopart
{
    protected:


        // pointer we will set to one of the thresholding functions
        typedef VectorXd (twopart::*thresh_func_2p_ptr)(VectorXd &value, VectorXd & penalty_factor, double &penalty, double &l1, double &denom);



        const MapMatd X, Xs;
        const MapVecd Z, S;
        VectorXd weights, weights_s, offset, offset_s, lambda;
        VectorXi groups, unique_groups;
        VectorXd group_weights;
        params P;

        double tau;
        int maxit, maxit_irls;
        double tol, tol_irls;

        bool intercept;
        std::string penalty;
        bool opposite_signs;

        bool need_to_compute_lambda;
        int nlambda;
        double lambda_min_ratio;

        int ngroups;
        std::vector<std::vector<int> > grp_idx;
        VectorXd eigenvals;


        int nobs, nobs_s, nvars;

        VectorXd beta, beta_old, beta_irls_old, beta_s, beta_s_old, beta_irls_s_old;
        MatrixXd beta_mat, beta_s_mat;

        VectorXd mu, mu_s, xbeta_cur, resid_cur, xbeta_s_cur, resid_s_cur, W, W_s;
        VectorXd deviance_vec, deviance_s_vec;



        double b0, b0_s, b0_old, b0_s_old, deviance, deviance_old, deviance_s, deviance_s_old, mult_1;

        VectorXd penalty_adjustment, ZZ;

        thresh_func_2p_ptr thresh_func;
        //VectorXd thresh_func(VectorXd &value, VectorXd & penalty_factor, double &penalty, double &l1, double &denom);


        VectorXd grad_func(const VectorXd &x_col,
                           const VectorXd &x_col_s,
                           VectorXd &xbeta,
                           VectorXd &xbeta_s);

        virtual void set_up_lambda();
        virtual void set_up_groups();

        // phi_j(v) function for cooperative lasso
        VectorXd phi_j_v(VectorXd & v, int & j);

        bool converged(const VectorXd& cur, const VectorXd& prev, const double& tolerance);

        bool converged_irls(double deviance, double deviance_prev, const double& tolerance);

        VectorXd compute_eigs_twopart();

        double soft_thresh(double & a, double & lambda);

        VectorXd block_soft_thresh_tp(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom);
        VectorXd coop_block_soft_thresh_tp(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom);




    public:
        twopart(const MapMatd & X_,
                const MapMatd & Xs_,
                const MapVecd & Z_,
                const MapVecd & S_,
                VectorXd & weights_,
                VectorXd & weights_s_,
                const MapVecd & offset_,
                const MapVecd & offset_s_,
                VectorXd & lambda_,
                const MapVeci & groups_,
                const MapVeci & unique_groups_,
                VectorXd & group_weights_,
                params & P_) :
        X(X_),
        Xs(Xs_),
        Z(Z_),
        S(S_),
        weights(weights_),
        weights_s(weights_s_),
        offset(offset_),
        offset_s(offset_s_),
        lambda(lambda_),
        groups(groups_),
        unique_groups(unique_groups_),
        group_weights(group_weights_),
        P(P_),
        tau(P.tau),
        maxit(P.maxit),
        maxit_irls(P.maxit_irls),
        tol(P.tol),
        tol_irls(P.tol_irls),
        intercept(P.intercept),
        penalty(P.penalty),
        opposite_signs(P.opposite_signs),
        need_to_compute_lambda(lambda.size() < 1),
        nlambda((need_to_compute_lambda) ? P.nlambda : lambda.size()),
        lambda_min_ratio(P.lambda_min_ratio),
        ngroups(unique_groups.size()),
        grp_idx(ngroups),
        eigenvals(ngroups),
        nobs(X.rows()),
        nobs_s(Xs.rows()),
        nvars(X.cols()),
        beta(VectorXd::Zero(nvars)),
        beta_old(VectorXd::Zero(nvars)),
        beta_irls_old(VectorXd::Zero(nvars)),
        beta_s(VectorXd::Zero(nvars)),
        beta_s_old(VectorXd::Zero(nvars)),
        beta_irls_s_old(VectorXd::Zero(nvars)),
        beta_mat(MatrixXd(nvars+1, nlambda)),
        beta_s_mat(MatrixXd(nvars+1, nlambda)),
        mu(nobs), mu_s(nobs_s), xbeta_cur(offset), resid_cur(nobs),
        xbeta_s_cur(offset_s), resid_s_cur(nobs_s), W(nobs), W_s(nobs_s),
        deviance_vec(nlambda), deviance_s_vec(nlambda)
        {}



        /*
        CD(const arma::mat& Xi, const arma::vec& yi, const Params& P);

        virtual double Objective(arma::vec & r, arma::sp_mat & B) = 0;

        virtual FitResult Fit() = 0;

        bool Converged();

        void SupportStabilized();

        static CD * make_CD(const arma::mat& Xi, const arma::vec& yi, const Params& P);
         */

        virtual void initialize();
        virtual VectorXi fit_path();

        virtual MatrixXd get_beta_z()  { return beta_mat; }
        virtual MatrixXd get_beta_s()  { return beta_s_mat; }
        virtual VectorXd get_lambda()  { return lambda; }
        virtual VectorXd get_dev_z()   { return deviance_vec; }
        virtual VectorXd get_dev_s()   { return deviance_s_vec; }
        virtual VectorXd get_pen_adj() { return penalty_adjustment; }
        virtual VectorXd get_eigs()    { return eigenvals; }


};

#endif
