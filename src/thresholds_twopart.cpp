

#include "thresholds.h"
#include "thresholds_twopart.h"

VectorXd block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &gamma, double &l2, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    double anorm = a.norm();
    double thresh_fact;

    for (int j = 0; j < alen; ++j)
    {
        thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / anorm) / denom;
        retval(j) = a(j) * thresh_fact;
    }

    return(retval);
}


VectorXd coop_block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &gamma, double &l2, double &denom)
{
    int alen = a.size();

    VectorXd retval(alen);

    double thresh_fact;

    for (int j = 0; j < alen; ++j)
    {
        VectorXd phi_j_vec = phi_j_v(a, j);
        thresh_fact = std::max(0.0, 1.0 - lambda * penalty_factor(j) / phi_j_vec.norm());
        retval(j) = a(j) * thresh_fact / denom;
    }

    return(retval);
}


thresh_func_twopart_ptr set_twopart_threshold_func(std::string & penalty)
{
    thresh_func_twopart_ptr thresh_func;

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
