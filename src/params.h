#ifndef PARAMS_H
#define PARAMS_H

#include "RcppEigen.h"

struct params
{

    double tau = 0.0;
    int maxit = 500;
    int maxit_irls = 500;
    double tol = 1e-5;
    double tol_irls = 1e-5;

    bool intercept = true;
    std::string penalty = "grp.lasso";
    bool opposite_signs = false;
    bool strongrule = true;

    int nlambda = 100;
    double lambda_min_ratio = 1e-3;
};

#endif
