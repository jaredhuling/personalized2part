
#' @export
mmbcd <- function(x, y,
                  groups = NULL,
                  group_weights = NULL,
                  weights = rep(1, NROW(x)),
                  penalty = c("grp.lasso", "coop.lasso"),
                  family = c("gaussian", "binomial", "gamma"),
                  algorithm = c("mm", "irls", "gaussian"),
                  lambda = NULL,
                  nlambda = 100L,
                  lambda_min_ratio = 1e-3,
                  maxit = 500,
                  tol = 1e-5,
                  maxit.irls = 25,
                  tol.irls = 1e-5,
                  intercept = TRUE)
{
    p <- NCOL(x)
    n <- NROW(x)

    penalty   <- match.arg(penalty)
    family    <- match.arg(family)
    algorithm <- match.arg(algorithm)

    if (is.null(groups))
    {
        groups <- 1:p
    }

    if (is.null(group_weights))
    {
        group_weights <- numeric(0)
    }
    if (is.null(lambda))
    {
        lambda <- numeric(0)
    }

    y <- setup_y(y, family)

    groups <- as.integer(groups)

    unique_groups    <- unique(groups)

    group_weights    <- as.double(group_weights)
    weights          <- as.double(weights)
    lambda           <- as.double(lambda)
    nlambda          <- as.integer(nlambda)
    maxit            <- as.integer(maxit)
    tol              <- as.double(tol)
    lambda_min_ratio <- as.double(lambda_min_ratio)
    intercept        <- as.logical(intercept)

    if (algorithm == "mm")
    {
        res <- mmbcd_cpp(X = x, Y = y, groups = groups,
                         unique_groups = unique_groups,
                         group_weights = group_weights,
                         weights = weights,
                         lambda = lambda, nlambda = nlambda,
                         lambda_min_ratio = lambda_min_ratio,
                         maxit = maxit, tol = tol, intercept = intercept,
                         penalty = penalty, family = family)
    } else if (algorithm == "irls")
    {
        res <- irls_mmbcd_cpp(X = x, Y = y, groups = groups,
                              unique_groups = unique_groups,
                              group_weights = group_weights,
                              weights = weights,
                              lambda = lambda, nlambda = nlambda,
                              lambda_min_ratio = lambda_min_ratio,
                              maxit = maxit, tol = tol,
                              maxit_irls = maxit.irls, tol_irls = tol.irls,
                              intercept = intercept,
                              penalty = penalty, family = family)
    } else
    {
        stopifnot(family == "gaussian")
        res <- mmbcd_gaussian_cpp(X = x, Y = y, groups = groups,
                                  unique_groups = unique_groups,
                                  group_weights = group_weights,
                                  weights = weights,
                                  lambda = lambda, nlambda = nlambda,
                                  lambda_min_ratio = lambda_min_ratio,
                                  maxit = maxit, tol = tol,
                                  intercept = intercept,
                                  penalty = penalty)
    }

    res
}


