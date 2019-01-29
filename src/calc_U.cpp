
#include "calc_U.h"

// function which sets the U function according to family distribution
U_func_ptr set_U_func(std::string & family)
{
    U_func_ptr U_func;

    if (family == "gaussian")
    {
        U_func = U_func_gaussian;
    } else if (family == "binomial")
    {
        U_func = U_func_binomial;
    } else if (family == "gamma")
    {
        U_func = U_func_gamma;
    }  else
    {
        U_func = U_func_gaussian;
    }

    return(U_func);
}

// function which sets the U function according to family distribution
U_intercept_func_ptr set_U_intercept_func(std::string & family)
{
    U_intercept_func_ptr U_int_func;

    if (family == "gaussian")
    {
        U_int_func = U_intercept_func_gaussian;
    } else if (family == "binomial")
    {
        U_int_func = U_intercept_func_binomial;
    } else if (family == "gamma")
    {
        U_int_func = U_intercept_func_gamma;
    } else
    {
        U_int_func = U_intercept_func_gaussian;
    }

    return(U_int_func);
}


/// GAUSSIAN

VectorXd U_func_gaussian(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    VectorXd U_vec = (x_subs.transpose() *
        (weights.array() * (y.array() - xbeta.array())).matrix()).array() / double(n);

    return(U_vec);
}

double U_intercept_func_gaussian(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = (weights.array() * (y.array() - xbeta.array())  ).matrix().sum() / double(n);
    return(U_val);
}


/// BINOMIAL

VectorXd U_func_binomial(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    VectorXd U_vec = (x_subs.transpose() *
        (weights.array() * (y.array() / (1.0 + (y.array() * xbeta.array()).exp()  ))).matrix()).array() / double(n);

    return(U_vec);
}

double U_intercept_func_binomial(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = ((weights.array() * (y.array() / (1.0 + (y.array() * xbeta.array()).exp()  ))).matrix()).sum() / double(n);

    return(U_val);
}

/// GAMMA

VectorXd U_func_gamma(MatrixXd &x_subs, const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    VectorXd U_vec = (x_subs.transpose() *
        (weights.array() * (y.array() * (-1.0 * xbeta.array()).exp() - 1.0  )).matrix()).array() / double(n);

    return(U_vec);
}

double U_intercept_func_gamma(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = ((weights.array() * (y.array() * (-1.0 * xbeta.array()).exp() - 1.0  )).matrix()).sum() / double(n);

    return(U_val);
}

