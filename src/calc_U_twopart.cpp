
#include "calc_U_twopart.h"
#include "calc_U.h"

// function which sets the U function according to family distribution
U_tp_func_ptr set_U_tp_func()
{
    U_tp_func_ptr U_func;

    return(U_func_twopart);
}

// function which sets the U function according to family distribution
U_tp_intercept_func_ptr set_U_tp_intercept_func(std::string & family)
{
    U_tp_intercept_func_ptr U_int_func;

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


/*
double U_intercept_func_gaussian(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = (weights.array() * (y.array() - xbeta.array())  ).matrix().sum() / double(n);
    return(U_val);
}
 */


/// TWO PART

VectorXd U_func_twopart(const VectorXd &x_col, const VectorXd &x_col_s,
                        const Eigen::Map<Eigen::VectorXd> &z,
                        const Eigen::Map<Eigen::VectorXd> &s,
                        VectorXd &weights,
                        VectorXd &weights_s,
                        VectorXd &xbeta,
                        VectorXd &xbeta_s,
                        int & n,
                        int & n_s)
{
    VectorXd U_vec(2);

    // binomial part
    U_vec(0) = (x_col.array() *
        (weights.array() * (z.array() / (1.0 + (z.array() * xbeta.array()).exp()  )))).matrix().sum() / double(n);

    // gamma part
    U_vec(1) = (x_col_s.array() *
        (weights_s.array() * (s.array() * (-1.0 * xbeta_s.array()).exp() - 1.0  ))).matrix().sum() / double(n_s);

    return(U_vec);
}

/*
double U_intercept_func_binomial(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = ((weights.array() * (y.array() / (1.0 + (y.array() * xbeta.array()).exp()  ))).matrix()).sum() / double(n);

    return(U_val);
}
 */

/*
double U_intercept_func_gamma(const Eigen::Map<Eigen::VectorXd> & y, VectorXd &weights, VectorXd &xbeta, int & n)
{
    double U_val = ((weights.array() * (y.array() * (-1.0 * xbeta.array()).exp() - 1.0  )).matrix()).sum() / double(n);

    return(U_val);
}
 */

