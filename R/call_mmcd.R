
#' @export
mmbcd <- function(x, y,
                  groups = NULL,
                  group_weights = NULL,
                  weights = rep(1, NROW(x)),
                  offset = NULL,
                  penalty = c("grp.lasso", "coop.lasso"),
                  family = c("gaussian", "binomial", "gamma"),
                  algorithm = c("mm", "irls", "gaussian"),
                  lambda = NULL,
                  nlambda = 100L,
                  lambda_min_ratio = ifelse(n < p, 0.05, 0.005),
                  tau = 0,
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

    if (is.null(weights))
    {
        weights <- rep(1, n)
    }

    if (is.null(offset))
    {
        offset <- rep(0, n)
    }

    stopifnot(length(groups) == p)
    stopifnot(length(weights) == n)
    stopifnot(length(offset) == n)

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

    if (length(group_weights) > 0)
    {
        stopifnot(length(group_weights) == length(unique_groups))
    }

    group_weights    <- as.double(group_weights)
    weights          <- as.double(weights)
    offset           <- as.double(offset)
    lambda           <- as.double(lambda)
    nlambda          <- as.integer(nlambda)
    maxit            <- as.integer(maxit)
    tol              <- as.double(tol)
    lambda_min_ratio <- as.double(lambda_min_ratio)
    intercept        <- as.logical(intercept)

    alpha <- 1
    alpha <- as.double(alpha[1])
    tau   <- as.double(tau[1])

    if (alpha < 0 | alpha > 1)
    {
        stop("alpha must be between 0 and 1")
    }

    if (tau < 0 | tau > 1)
    {
        stop("tau must be between 0 and 1")
    }

    if (algorithm == "mm")
    {
        res <- mmbcd_cpp(X = x, Y = y, groups = groups,
                         unique_groups = unique_groups,
                         group_weights = group_weights,
                         weights = weights, offset = offset,
                         lambda = lambda, nlambda = nlambda,
                         lambda_min_ratio = lambda_min_ratio,
                         alpha = alpha, tau = tau,
                         maxit = maxit, tol = tol, intercept = intercept,
                         penalty = penalty, family = family)
    } else if (algorithm == "irls")
    {
        res <- irls_mmbcd_cpp(X = x, Y = y, groups = groups,
                              unique_groups = unique_groups,
                              group_weights = group_weights,
                              weights = weights, offset = offset,
                              lambda = lambda, nlambda = nlambda,
                              lambda_min_ratio = lambda_min_ratio,
                              alpha = alpha, tau = tau,
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


