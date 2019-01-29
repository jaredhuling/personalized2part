

#include "thresholds.h"

VectorXd block_soft_thresh(VectorXd & a, double & lambda, double &gamma, double &l2, double &denom)
{
    double thresh_fact = std::max(0.0, 1.0 - lambda / a.norm());
    VectorXd retval = a.array() * thresh_fact / denom;

    return(retval);
}

VectorXd phi_j_v(VectorXd & v, int & j)
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

VectorXd coop_block_soft_thresh(VectorXd & a, double & lambda, double &gamma, double &l2, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    double thresh_fact;

    for (int j = 0; j < alen; ++j)
    {
        VectorXd phi_j_vec = phi_j_v(a, j);
        thresh_fact = std::max(0.0, 1.0 - lambda / phi_j_vec.norm());
        retval(j) = a(j) * thresh_fact / denom;
    }

    return(retval);
}


double soft_thresh(double & a, double & lambda)
{
    double thresh_fact = std::max(0.0, 1.0 - lambda / std::abs(a) );
    double retval = a * thresh_fact;

    return(retval);
}

thresh_func_ptr set_threshold_func(std::string & penalty)
{
    thresh_func_ptr thresh_func;

    if (penalty == "grp.lasso")
    {
        thresh_func = block_soft_thresh;
    } else if (penalty == "coop.lasso")
    {
        thresh_func = coop_block_soft_thresh;
    } else
    {
        thresh_func = block_soft_thresh;
    }

    return(thresh_func);
}
