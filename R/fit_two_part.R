

#' @export
mmbcd_twopart <- function(x, z,
                          x.s, s,
                          penalty.factor   = NULL,
                          weights          = rep(1, NROW(x)),
                          weights.s        = rep(1, NROW(x.s)),
                          penalty          = c("grp.lasso", "coop.lasso"),
                          algorithm        = c("mm", "irls"),
                          lambda           = NULL,
                          nlambda          = 100L,
                          lambda_min_ratio = 1e-3,
                          intercept        = TRUE,
                          maxit            = 500,
                          tol              = 1e-5,
                          maxit.irls       = 25,
                          tol.irls         = 1e-5)
{
    p <- NCOL(x)
    n <- NROW(x)

    algorithm <- match.arg(algorithm)
    penalty   <- match.arg(penalty)

    if (is.null(penalty.factor))
    {
        penalty.factor <- numeric(0)
    }
    if (is.null(lambda))
    {
        lambda <- numeric(0)
    }


    z <- setup_y(z, "binomial")
    s <- setup_y(s, "gamma")

    groups           <- as.integer(rep(1:NCOL(x), 2))
    unique_groups    <- unique(groups)

    penalty.factor   <- as.double(penalty.factor)
    weights          <- as.double(weights)
    weights.s        <- as.double(weights.s)
    lambda           <- as.double(lambda)
    nlambda          <- as.integer(nlambda)
    maxit            <- as.integer(maxit)
    tol              <- as.double(tol)
    maxit.irls       <- as.integer(maxit.irls)
    tol.irls         <- as.double(tol.irls)
    lambda_min_ratio <- as.double(lambda_min_ratio)
    intercept        <- as.logical(intercept)

    if (algorithm == "mm")
    {
        res <- mmbcd_twopart_cpp(X = x, Z = z,
                                 Xs = x.s, S = s,
                                 groups = groups,
                                 unique_groups = unique_groups,
                                 group_weights = penalty.factor,
                                 weights = weights,
                                 weights_s = weights.s,
                                 lambda = lambda, nlambda = nlambda,
                                 lambda_min_ratio = lambda_min_ratio,
                                 maxit = maxit, tol = tol, intercept = intercept,
                                 penalty = penalty)
    } else
    {
        res <- irls_mmbcd_twopart_cpp(X = x, Z = z,
                                      Xs = x.s, S = s,
                                      groups = groups,
                                      unique_groups = unique_groups,
                                      group_weights = penalty.factor,
                                      weights = weights,
                                      weights_s = weights.s,
                                      lambda = lambda, nlambda = nlambda,
                                      lambda_min_ratio = lambda_min_ratio,
                                      maxit = maxit, tol = tol,
                                      maxit_irls = maxit.irls, tol_irls = tol.irls,
                                      intercept = intercept,
                                      penalty = penalty)
    }

    res
}
