

#include "thresholds.h"
#include "thresholds_twopart.h"

VectorXd block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
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


VectorXd coop_block_soft_thresh(VectorXd & a, VectorXd & penalty_factor, double & lambda, double &l1, double &denom)
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
